# functions_rcpp_maps.R
# Main functions for running the DynaMICE model based on Rcpp files
# update: 2022/11/22

# ------------------------------------------------------------------------------
# Measles vaccine impact model
#   To estimate the health impact of measles vaccination for a given set of
#   countries.
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
## Execute the Rcpp measles model for a country
#-------------------------------------------------------------------------------
# Notes for change
# 1. fix spin-up period = 100000
# 2. change inputs for SIAs based on delivery strategies
#-------------------------------------------------------------------------------
runCountry_rcpp <- function (
  #variables specific for loop
  iso3,
  years,

  #infection dynamic variables
  vaccination,
  using_sia,
  parms_rcpp,

  #input data
  c_coverage_routine,
  c_coverage_sia,
  c_timeliness,
  contact_mat,
  c_contact,
  c_rnought,
  c_population,

  # dynaMICE model options
  tstep,
  save_scenario,
  foldername,

  r

) {

  # adjusted filename
  if (iso3 == "XK"){iso3 <- "XKX"}

  # country specific timeliness curve
  country_timeliness <- c_timeliness [!is.na(age), timeliness]
  timeliness_ages    <- c_timeliness [!is.na(age), age]

  # expand 0-2 years old to weekly age strata
  s         <- 52 # number of finer stages within an age band (weekly ages, so 52)
  jt        <- 3  # how many ages to expand to s (or to weeks)
  beta_full <- matrix (0, ncol = 254, nrow = 254)

  beta_full[(1:(s*jt)), (1:(s*jt))] <- expandMatrix (
    A = c_contact [1:jt, 1:jt]/s,  # needs to be divided by 52 so that the mean total number of contacts stays the same
    expand_rows =  s, expand_cols =  s,
    rescale_rows = FALSE, rescale_cols = FALSE)

  beta_full[1:(s*jt),((s*jt)+1):(ncol(beta_full))] <- expandMatrix(
    A = c_contact [1:jt,(jt+1):ncol(c_contact)],
    expand_rows = s, expand_cols = 1,
    rescale_rows = F, rescale_cols = F)

  beta_full[((s*jt)+1):(nrow(beta_full)), 1:(s*jt)] <- expandMatrix(
    A = c_contact [(jt+1):nrow(c_contact),1:jt]/s,  # adjust to ensure the mean total number of contacts stays the same
    expand_rows = 1, expand_cols = s,
    rescale_rows = F, rescale_cols = F)

  beta_full[((s*jt)+1):(nrow(beta_full)), ((s*jt)+1):(ncol(beta_full))] <-
    c_contact [(jt+1):nrow(c_contact),(jt+1):ncol(c_contact)]

  beta_full_unadj <- beta_full

  # infection rate under target R0
  beta_tstep     <- c_rnought*parms_rcpp$gamma

  # setup inputs for DynaMICE Rcpp model
  n_years        <- length(years)

  #y_out          <- array(0, c(254, 14, n_years))   # numbers of age groups, compartments, years
  case_out       <- array(0, c(254, n_years))
  case0d_out     <- array(0, c(254, n_years))
  pop_out        <- array(0, c(254, n_years))
  popS_out       <- array(0, c(254, n_years))
  popS0d_out     <- array(0, c(254, n_years))
  popS1d_out     <- array(0, c(254, n_years))
  popR_out       <- array(0, c(254, n_years))
  doseRI1_out    <- array(0, c(254, n_years))
  doseRI2_out    <- array(0, c(254, n_years))
  doseSIAc_out   <- array(0, c(254, n_years))
  doseSIAde1_out <- array(0, c(254, n_years))
  doseSIAde2_out <- array(0, c(254, n_years))
  doseSIAf_out   <- array(0, c(254, n_years))
  imm1pop_out    <- array(0, c(254, n_years))
  imm2pop_out    <- array(0, c(254, n_years))

  init_Comp     <- matrix(0, 254, 14)
  init_Comp[,2] <- 0.95
  init_Comp[,3] <- 0.05

  spinup_yr <- 100  # duration in years of spinup period

  writelog ("log", paste0(iso3, "; Started model run"))
  # run model by yearly input
  for (y in years) {

    #old script groups those aged 70-80, but division is by actual population size
    pop.vector <- c_population[year == y, value]

    # first expand polymod matrix (contact_tstep) and population vector and
    # then divide by population sizes, otherwise it doesn't work.
    pop.vector_full <- c(rep(pop.vector[1:jt]/s, each = s), pop.vector[(jt +1):length(pop.vector)])

    # change zero values to 1 to avoid division by 0
    pop.vector_full[pop.vector_full==0] <- 1

    # generate contact matrices for DynaMICE using uniform age structure
    # adjust contact reciprocity by population structure of each year
    if (contact_mat == "unimix") {

      beta_full <- beta_tstep * beta_full_unadj

      } else if (contact_mat == "prpmix") {

        beta_full <- beta_tstep * matrix (rep(pop.vector_full/sum(pop.vector_full), each = 254), nrow = 254)

      } else {

        beta_full <- matrix(0, 254, 254)
        beta_full_R0 <- matrix(0, 254, 254)
        for (i in 1:254){
          for (j in 1:254) {
            beta_full [i, j] <- (beta_full_unadj[i, j] * pop.vector_full[i] +
                                   beta_full_unadj[j, i] * pop.vector_full[j])/(2*pop.vector_full[i])

            # calculate infection rate of the contact matrix based on the R0 definition
            # transform from "contactees per contactor" to "contactors per contactee"
            # This step can be skipped, since the largest eigenvalue remains unchanged.
            beta_full_R0 [i, j] <-  beta_full [i, j] * (pop.vector_full[i] / pop.vector_full[j])
          }
        }
        # make sure the contact matrix to represent target R0
        beta_full <- (beta_tstep / Re(eigen(beta_full_R0, only.values=T)$values[1])) * beta_full
      }

      # run spin-up period
      if (y == years[1]){
        outp <- rcpp_spinup_maps (init_Comp, parms_rcpp, beta_full, pop.vector_full, spinup_yr)
        out_Comp <- outp$out_Comp

      if (vaccination >= 1) {

        # Maximum coverage can (obviously) only be 100%
        # To estimate proportion that is vaccinated at each week, we first calculate the number of individuals remaining susceptible
        # Then we calculate the number of individuals that should be vaccinated each week, in order to remain 1 - coverage susceptibles at the end of the timeliness data
        # In essence, this becomes the inverse of the cumulative timeliness curve
        cycov <- c_coverage_routine [year == y & vaccine == "MCV1", coverage]

        # Not use the following adjustment for coverage, as it leads to reduced coverage estimates
        # cycov <- c_coverage_routine [year == y & vaccine == "MCV1", coverage] /
        #   c_timeliness [is.na(age), prop_final_cov]

        if (length(cycov) == 0) {   # check if vaccine is not yet introduced and thereby, coverage value missing for this year
          cycov <- 0
        } else if (is.na(cycov)) {
          cycov <- 0
        }

        country_year_timeliness_mcv1 <- 1 - min(cycov,1) * country_timeliness

        country_year_timeliness_mcv1 <- -diff(country_year_timeliness_mcv1) /
          (country_year_timeliness_mcv1 [1:(length(country_year_timeliness_mcv1)-1)])

        # Timeliness is reported by week in the first year, and by month in the second year. Assume there is no vaccination in between
        country_year_timeliness_mcv1 [is.na(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1 [is.nan(country_year_timeliness_mcv1)] <- 0
        country_year_timeliness_mcv1_allages <- rep(0, 254)
        country_year_timeliness_mcv1_allages [round(timeliness_ages)] <- country_year_timeliness_mcv1

      } else {
        country_year_timeliness_mcv1_allages <- rep(0, 254)
      }

      if(vaccination == 2){
        country_year_mcv2 <- c_coverage_routine [year == y & vaccine == "MCV2", coverage]
      } else {
        country_year_mcv2 <- 0
      }

      if ( is.na(country_year_mcv2) || length (country_year_mcv2) == 0 ) {
        country_year_mcv2 <- 0
      }

      # SIA inputs
      if ((using_sia == 1) & (nrow (c_coverage_sia) > 0)) {

        setorder (c_coverage_sia, country_code, year, strategy)
        sia_targetpop <- tolower (stringr::str_sub (c_coverage_sia [year == y, strategy], 1, 1))
        sia_targetpop <- which (letters %in% sia_targetpop) # convert to numbers
        # sia_targetpop <- ifelse (sia_targetpop == "C" | sia_targetpop == "F", 0,
        #                          ifelse (sia_targetpop == "D", 1, 2))
        sia_input <- list (
          sia_rounds = length(sia_targetpop),
          a0 = c_coverage_sia [year == y, a0],
          a1 = c_coverage_sia [year == y, a1],
          siacov = c_coverage_sia [year == y, coverage],
          poptype = as.integer(sia_targetpop)
        )
      } else {
        sia_input <- list (
          sia_rounds = 0,
          a0 = integer(0),
          a1 = integer(0),
          siacov = numeric(0),
          poptype = integer(0)
        )
      }

      # run main period
      t_start <- spinup_yr*tstep + (y-years[1])*tstep + 1
      outp <- rcpp_vaccine_oney_maps (out_Comp,
                                      parms_rcpp,
                                      sia_input,
                                      beta_full,
                                      pop.vector_full,
                                      country_year_timeliness_mcv1_allages,
                                      country_year_mcv2,
                                      t_start)
      out_Comp <- outp$out_Comp


      if (length (which (out_Comp < 0)) > 0) {
        print (paste0 ("age group: ", which (out_Comp < 0, arr.ind = T)[,1]))
        print (paste0 ("state: ", which (out_Comp < 0, arr.ind = T)[,2]))
               }

      case_out       [, (y-years[1])+1] <- outp$cases*pop.vector_full        # new infections adding to I, V1I, V2I, V3I
      case0d_out     [, (y-years[1])+1] <- outp$cases0d*pop.vector_full      # new infections adding to I
      pop_out        [, (y-years[1])+1] <- rowSums(out_Comp[, 1:13])*pop.vector_full  # all compartments
      popS_out       [, (y-years[1])+1] <- rowSums(out_Comp[, c(2, 5, 8, 11)])*pop.vector_full
      popS0d_out     [, (y-years[1])+1] <- out_Comp[, 2]*pop.vector_full
      popS1d_out     [, (y-years[1])+1] <- out_Comp[, 5]*pop.vector_full
      popR_out       [, (y-years[1])+1] <- rowSums(out_Comp[, c(4, 7, 10, 13)])*pop.vector_full
      doseRI1_out    [, (y-years[1])+1] <- outp$doseRI1*pop.vector_full
      doseRI2_out    [, (y-years[1])+1] <- outp$doseRI2*pop.vector_full
      doseSIAc_out   [, (y-years[1])+1] <- outp$doseSIAc*pop.vector_full
      doseSIAde1_out [, (y-years[1])+1] <- outp$doseSIAde1*pop.vector_full
      doseSIAde2_out [, (y-years[1])+1] <- outp$doseSIAde2*pop.vector_full
      doseSIAf_out   [, (y-years[1])+1] <- outp$doseSIAf*pop.vector_full
      imm1pop_out    [, (y-years[1])+1] <- rowSums(out_Comp[, 5:7])*pop.vector_full
      imm2pop_out    [, (y-years[1])+1] <- rowSums(out_Comp[, 8:13])*pop.vector_full

      #if(y %% 20 == 0) {print (paste0 ('year ', y, ' finished'))}
  }

    writelog ("log", paste0 (iso3,"; Finished model run"))

    saveRDS (list(cases      = rbind (colSums(      case_out[1:52,]), colSums(      case_out[53:104,]), colSums(      case_out[105:156,]),       case_out[157:254,]),
                  cases0d    = rbind (colSums(    case0d_out[1:52,]), colSums(    case0d_out[53:104,]), colSums(    case0d_out[105:156,]),     case0d_out[157:254,]),
                  pops       = rbind (colSums(       pop_out[1:52,]), colSums(       pop_out[53:104,]), colSums(       pop_out[105:156,]),        pop_out[157:254,]),
                  popSus     = rbind (colSums(      popS_out[1:52,]), colSums(      popS_out[53:104,]), colSums(      popS_out[105:156,]),       popS_out[157:254,]),
                  popSus0d   = rbind (colSums(    popS0d_out[1:52,]), colSums(    popS0d_out[53:104,]), colSums(    popS0d_out[105:156,]),     popS0d_out[157:254,]),
                  popSus1d   = rbind (colSums(    popS1d_out[1:52,]), colSums(    popS1d_out[53:104,]), colSums(    popS1d_out[105:156,]),     popS1d_out[157:254,]),
                  popRec     = rbind (colSums(      popR_out[1:52,]), colSums(      popR_out[53:104,]), colSums(      popR_out[105:156,]),       popR_out[157:254,]),
                  doseRI1    = rbind (colSums(   doseRI1_out[1:52,]), colSums(   doseRI1_out[53:104,]), colSums(   doseRI1_out[105:156,]),    doseRI1_out[157:254,]),
                  doseRI2    = rbind (colSums(   doseRI2_out[1:52,]), colSums(   doseRI2_out[53:104,]), colSums(   doseRI2_out[105:156,]),    doseRI2_out[157:254,]),
                  doseSIAc   = rbind (colSums(  doseSIAc_out[1:52,]), colSums(  doseSIAc_out[53:104,]), colSums(  doseSIAc_out[105:156,]),   doseSIAc_out[157:254,]),
                  doseSIAde1 = rbind (colSums(doseSIAde1_out[1:52,]), colSums(doseSIAde1_out[53:104,]), colSums(doseSIAde1_out[105:156,]), doseSIAde1_out[157:254,]),
                  doseSIAde2 = rbind (colSums(doseSIAde2_out[1:52,]), colSums(doseSIAde2_out[53:104,]), colSums(doseSIAde2_out[105:156,]), doseSIAde2_out[157:254,]),
                  doseSIAf   = rbind (colSums(  doseSIAf_out[1:52,]), colSums(  doseSIAf_out[53:104,]), colSums(  doseSIAf_out[105:156,]),   doseSIAf_out[157:254,]),
                  imm1pop    = rbind (colSums(   imm1pop_out[1:52,]), colSums(   imm1pop_out[53:104,]), colSums(   imm1pop_out[105:156,]),    imm1pop_out[157:254,]),
                  imm2pop    = rbind (colSums(   imm2pop_out[1:52,]), colSums(   imm2pop_out[53:104,]), colSums(   imm2pop_out[105:156,]),    imm2pop_out[157:254,])),
             file = paste0 ("outcome/", save_scenario, "/", foldername, "/", iso3, ".RDS")
             )
  }
}
# end of function -- runCountry_rcpp
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
## Run the Rcpp measles model for a selected vaccination strategy
#-------------------------------------------------------------------------------
# Notes for change
# 1. use India data on rnought and timeliness for all countries
#-------------------------------------------------------------------------------
runScenario_rcpp <- function (vaccine_coverage_folder    = "",
                              scenario_name,
                              save_scenario,
                              burden_estimate_folder,
                              countries                  = "all",
                              vaccination,  # Whether children are vaccinated. 0: No vaccination; 1: Only MCV1; 2: MCV1 and MCV2
                              using_sia,    # Whether supplementary immunization campaigns are used. 0: no SIA; 1: with SIA
                              contact_mat,                             # contact matrix: "prpmix","polymod","syn"
                              sim_years                                # calendar years for simulation
) {

  ages           <- c (0:100)   # Numeric vector with age-strata that are reported (Model ALWAYS models 101 age-groups; 0-2yo in weekly age-strata, and 3-100yo in annual age-strata)

  # --------------------------------------------------------------------------
  # input data file names
  # --------------------------------------------------------------------------
  # vaccine coverage file for routine immunisation
  data_coverage_routine <- paste0 (vaccine_coverage_folder,
                                   "routine_",
                                   scenario_name,
                                   ".csv")

  # separate vaccine coverage file for SIA (supplementary immunisation activities)
  data_coverage_sia <- paste0 (vaccine_coverage_folder,
                               "sia_",
                               scenario_name,
                               ".csv")


  # --------------------------------------------------------------------------
  # Advanced options
  # --------------------------------------------------------------------------
  dinf		    <- 14				                   # duration of infection (days)
  tstep			  <- 1000				                 # number of timesteps in a year

  # age-dependent vaccine efficacy for first dose, based on a linear model (Hughes et al. 2020)
  ve1_intcp   <- 0.64598                     # intercept of the linear model
  ve1_slope   <- 0.01485                     # slope of the linear model, per month of age

  ve2plus     <- 0.98                        # vaccine efficacy for two and moer doses

  age_ve1 <- rep (NA, 254)
  age_ve1[1:(3*52)]      <- ve1_intcp + ve1_slope*((1:(3*52))/52)*12
  age_ve1[(3*52+1):254]  <- ve1_intcp + ve1_slope*(4:101)*12
  age_ve1 <- ifelse (age_ve1 >= ve2plus, ve2plus, age_ve1)

  # parameters to include for Rcpp functions
  parms_rcpp <- list (gamma         = 1 / (dinf * tstep/365),    # rate of losing infectivity
                      tstep         = tstep,
                      amp           = 0.05,		                   # amplitude for seasonality
                      ve1           = age_ve1,
                      ve2plus       = ve2plus)

  # create folders with correct name for in- and output data if not yet exists
  # typically foldername should not exist - but may come in handy when only processing results
  if ( !exists("foldername_analysis") ) {

    foldername <- paste0 (
      format(Sys.time(),format="%Y%m%d"),
      "_v",
      vaccination,
      "_s",
      using_sia,
      "_deter"
    )

    dir.create(
      file.path(
        paste0(
          getwd(),
          "/outcome/", save_scenario, "/",
          foldername
        )
      ), recursive = T
    )

  } else {
    foldername <- foldername_analysis
  }

  # log
  writelog ("log", paste0("Main; ", scenario_name, " started"))


  # --------------------------------------------------------------------------
  # prepare input data
  # --------------------------------------------------------------------------
  # read_data
  coverage_routine	<- copy (fread (data_coverage_routine))[country_code %in% countries & year %in% sim_years]
  coverage_sia		  <- copy (fread (data_coverage_sia))[country_code %in% countries & year %in% sim_years]
  timeliness  		  <- setDT (data_timeliness)[country_code %in% countries]
  rnought	    		  <- setDT (data_r0)[country_code %in% countries]
  population  		  <- setDT (data_pop)[country_code %in% countries]

  # Process contact matrices
  contact_list <- switch (contact_mat,
                          "syn"     = sapply (countries,
                                              function(cty){data_contact_syn[[cty]]},
                                              simplify = FALSE, USE.NAMES = TRUE)
  )


  # ----------------------------------------------------------------------------
  # Run model
  # ----------------------------------------------------------------------------
  for (ii in 1:length(countries)){
    iso3    <- countries[ii]
    out_run <- runCountry_rcpp (iso3               = iso3,
                                years              = as.numeric (sim_years),
                                vaccination        = vaccination,
                                using_sia          = using_sia,
                                parms_rcpp         = parms_rcpp,
                                c_coverage_routine = coverage_routine[country_code == iso3,],
                                c_coverage_sia     = coverage_sia[country_code == iso3 & coverage != 0,],
                                c_timeliness       = timeliness[country_code == iso3,],
                                contact_mat        = contact_mat,
                                c_contact          = contact_list[[iso3]],
                                c_rnought          = rnought[country_code == iso3, r0],
                                c_population       = population[country_code == iso3,],
                                tstep              = tstep,
                                save_scenario      = save_scenario,
                                foldername         = foldername
    )
  }

  writelog ("log", paste0 ("Main; finished"))
  return ()

} # end of function -- runScenario_rcpp
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
## Output a csv file for burden estimates
# ------------------------------------------------------------------------------
output_burden_estimate <- function (save_scenario,
                                    foldername,
                                    sim_years,
                                    cfr_option,
                                    burden_estimate_folder,
                                    burden_estimate_file,
                                    vimc_scenario,
                                    portnoy_scenario
) {

  # ----------------------------------------------------------------------------
  # merge and process results
  # ----------------------------------------------------------------------------
  ages     <- c(0:100)
  years    <- as.numeric (sim_years)

  # merge results
  output_files <- list.files (path = paste0 ("outcome/", save_scenario, "/", foldername, "/"),
                              recursive = T, full.names = T)
  all_runs  <- rbindlist (lapply (output_files, function (filename, ...) {
    res <- withCallingHandlers (
      readRDS (filename),
      warning = function(w) {warning(w, filename); }
    )
    filefinal <- stringr::str_extract (filename, "[^/]+$")      # remove path but keep filename
    res2 <- data.table (country     = rep (stringr::str_sub(filefinal, 1, 3), 101*length(years)),
                        year        = rep (years, each = 101),
                        age         = rep (0:100, length(years)),
                        cases       = as.vector(res$cases),
                        cases0d     = as.vector(res$cases0d),
                        cohort_size = as.vector(res$pops),
                        popSus      = as.vector(res$popSus),
                        popSus0d    = as.vector(res$popSus0d),
                        popSus1d    = as.vector(res$popSus1d),
                        popRec      = as.vector(res$popRec),
                        doseRI1     = as.vector(res$doseRI1),
                        doseRI2     = as.vector(res$doseRI2),
                        doseSIAc    = as.vector(res$doseSIAc),
                        doseSIAde1  = as.vector(res$doseSIAde1),
                        doseSIAde2  = as.vector(res$doseSIAde2),
                        doseSIAf    = as.vector(res$doseSIAf),
                        imm1pop     = as.vector(res$imm1pop),
                        imm2pop     = as.vector(res$imm2pop)
    )
    return (res2)
  }))

  # ----------------------------------------------------------------------------
  # add data column for remaining life expectancy
  # ----------------------------------------------------------------------------
  lexp_remain <- tailor_data_lexp_remain(unique(all_runs$country), years)

  cols_mea <- c("cases", "cases0d", "cohort_size",
                "popSus", "popSus0d", "popSus1d", "popRec",
                "doseRI1", "doseRI2", "doseSIAc",
                "doseSIAde1", "doseSIAde2", "doseSIAf",
                "imm1pop", "imm2pop")

  save.cols <- c("i.country", "year", "age", cols_mea, "value")

  all_runs <- lexp_remain [all_runs, ..save.cols,
                           on = .(country_code = country,
                                  age          = age,
                                  year         = year) ]

  # rename column names for output
  setnames (all_runs,
            old = c("i.country", "value"      ),
            new = c("country"  , "remain_lexp"))

  # ----------------------------------------------------------------------------
  # estimate deaths using Portnoy's method
  # ----------------------------------------------------------------------------
  # read CFRs (rates between 0 and 1)
  cfr = setDT (data_cfr_portnoy)

  # cfr estimates are for 2000 to 2030
  # if cfr estimates are required for years below or above this range, then
  # for years below 2000, set cfr estimates of year 2000
  # for years above 2030, set cfr estimates of year 2030

  # find minimum and maximum year
  min_year = min (all_runs [, year])
  max_year = max (all_runs [, year])

  # set cfrs for years before 2000 to cfr values of year 2000
  if (min_year < 2000) {
    cfr.year <- rbindlist (lapply (min_year:1999, function(i) copy (cfr [year == 2000, ])[, year := i]))
    cfr      <- rbind     (cfr.year, cfr, use.names = TRUE)
  }

  # set cfrs after 2030 to 2030
  if (max_year > 2030) {
    cfr.year <- rbindlist (lapply (2031:max_year, function(i) copy (cfr [year == 2030, ])[, year := i]))
    cfr      <- rbind     (cfr, cfr.year, use.names = TRUE)
  }

  # rename columns -- cfr of vimc_scenario and portnoy_scenario
  # cfrs for under 5 (< 5) and over 5 (>= 5) years
  setnames (x = cfr,
            old = c(paste0 ("cfr_under5_", vimc_scenario, "_", portnoy_scenario),
                    paste0 ("cfr_over5_" , vimc_scenario, "_", portnoy_scenario)),
            new = c("under5", "over5"))

  # add CFR data column to burden estimates
  save.cols <- c("year", "age", "country", cols_mea, "over5", "under5", "remain_lexp")

  burden <- cfr [all_runs, ..save.cols,
                 on = .(country_code = country, year = year) ]

  # estimate deaths for ages under 5 years
  burden [age < 5, deaths := cases * under5]

  # estimate deaths for ages over 5 years
  burden [age >= 5, deaths := cases * over5]


  # ----------------------------------------------------------------------------
  # DALYs
  # ----------------------------------------------------------------------------
  # calculate dalys = (ylds) + (ylls)
  burden [, dalys := ((cases - deaths) * 0.002) + (deaths * remain_lexp)]

  # adjust columns for output
  save.cols <- c("year", "age", "country",
                 cols_mea [c(2,1)], "deaths", "dalys", cols_mea [3:length(cols_mea)])
  burden <- subset (burden, select = save.cols)

  # append/suffix cfr_option to the end of filename
  updated_burden_estimate_file <- stringr::str_replace (string      = burden_estimate_file,
                                                        pattern     = ".csv",
                                                        replacement = paste0 ("_", cfr_option, ".csv")
  )

  # save updated burden estimate file (cases + deaths) to file
  # cfr_option is also the name of the subfolder
  fwrite (x    = burden,
          file = paste0 (burden_estimate_folder, "/",
                         updated_burden_estimate_file))
  return ()

} # end of function -- estimate_deaths_dalys
# ------------------------------------------------------------------------------

