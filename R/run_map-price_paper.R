# run_map-price_paper.R
# find country-specific thresholds for MR-MAP price
# modified from output_results_paper.R
# update: 2022/11/29

library(data.table)
library(stringr)
library(tidyr)

rm (list = ls())
# # load input information:
# # source("R/process_country_input_paper.R")
# # (1) dat_input: MAP intro year, unit costs, coverage, reported cases, funding groups, WHO region, Gavi-eligibility
# # (2) dat_oppcost_income: income-level health opportunity costs
# # (3) dat_oppcost_ctry:  country-level health opportunity costs
# # (4) data_pop: UNWPP 2019 population size by age and year
# # (5) dat_income: country list by income group

load (file = "inputs/dat_input.rds")
load (file = "inputs/dat_oppcost_income.rds")
load (file = "inputs/dat_oppcost_ctry.rds")
load (file = "data/data_pop_maps.rda")
load (file = "inputs/dat_income.rds")


# ------------------------------------------------------------------------------
## load coverage and result files by scenarios
# ------------------------------------------------------------------------------
# list of scenarios
scns <- c("high-base",  # (1) without MR-MAPs (High)
          "high-maps",  # (2) with MR-MAPs (High)
          "high-accl",  # (3) accelerated MR-MAP intro (High)
          "low-base",   # (4) without MR-MAPs (Low)
          "low-maps",   # (5) with MR-MAPs (Low)
          "low-accl")   # (6) accelerated MR-MAP intro (Low)
income_names <- c("Low income", "Lower middle income", "Upper middle income")

# exclude countries with low measles burden
rm_ctries <- c("ALB", "BWA", "CPV", "DJI", "ERI", "FJI", "GEO", "GIN", "GNB", "JOR",
               "KAZ", "KGZ", "MAR", "MUS", "MYS", "RWA", "SSD", "THA", "TON", "VUT")
anl_countries <- dat_income [!(income_g == "High income" | country_code %in% rm_ctries),
                             country_code] # 70

# number of cases between 2030-2040
file_case  <- NULL
for (scname in scns){
  scn_case <- fread (paste0 ("burden_estimate/burden_", scname, "_Portnoy.csv"))
  scn_case [, comp := scname]
  file_case <- rbind (file_case,
                      scn_case [country %in% anl_countries & year %in% 2030:2040])
}
annual_burden <- file_case [, lapply(.SD, sum),
                            .SDcols = cases:doseSIAf,
                            by = c("country", "year", "comp")]

# get MCV1 & MCV2 coverages and their incremental changes compared to 2020
file_cov <- NULL
sel_cols <- c("vaccine", "strategy", "country_code", "year", "coverage")
for (scname in scns){
  input_scn_cov <- rbind (fread (paste0 ("vac_coverage_maps/routine_", scname, ".csv")) [, ..sel_cols],
                          fread (paste0 ("vac_coverage_maps/sia_", scname, ".csv")) [, ..sel_cols])
  input_scn_cov <- copy (input_scn_cov [country_code %in% anl_countries]) [, scenario := scname]
  for (stgname in c("A","B")){
    # set up baseline
    scn_cov <- input_scn_cov [str_sub(strategy,1,1) == stgname & year %in% 2030:2040]
    scn_cov_base <- input_scn_cov [str_sub(strategy,1,1) == stgname & year == 2020,
                                   c("country_code", "coverage", "scenario")]
    scn_cov [scn_cov_base, on = .(country_code, scenario), cov2020 := i.coverage]

    # add in coverage for HTR and MOV populations (strategy pairs: A-D(MCV1), B-E(MCV2))
    scn_cov_add <- input_scn_cov [str_sub(strategy,1,1) == ifelse (stgname == "A", "D", "E") &
                                    year %in% 2030:2040,
                                  c("country_code", "coverage", "scenario", "year")]
    scn_cov [scn_cov_add, on = .(country_code, scenario, year), covadd := i.coverage]
    scn_cov [is.na(covadd), covadd := 0]
    scn_cov [, covall := ifelse (coverage + covadd > 1, 1, coverage + covadd)]

    file_cov <- rbind (file_cov, scn_cov)
  }
}


# ------------------------------------------------------------------------------
## calculate coverage- and burden-dependent costs
# ------------------------------------------------------------------------------
# in this subset, strategy "A" = vaccine "MCV1", and strategy "B" = vaccine "MCV2"
file_cov <- setDT (pivot_wider (copy (file_cov [, .(scenario, vaccine, country_code,
                                                    year, covall, cov2020)]),
                                names_from = vaccine,
                                names_sep = "_",
                                values_from = covall:cov2020))
file_cov [, `:=` (inc_covMCV1 = covall_MCV1 - cov2020_MCV1,
                  inc_covMCV2 = covall_MCV2 - cov2020_MCV2)]
file_cov [inc_covMCV1 < 0, inc_covMCV1 := 0] # no cost reduction when coverage is less than the baseline
file_cov [inc_covMCV2 < 0, inc_covMCV2 := 0]

# merge coverage input data
annual_burden <-  annual_burden [file_cov,
                                 on = .(country = country_code,
                                        year = year,
                                        comp = scenario)]
dat_anal <- annual_burden [dat_input [country_code %in% anl_countries],
                           on = .(country = country_code, comp = comp)]
setnames (x = dat_anal, old = "i.country", new = "country_name")

# adjust market penetration rate based on scenario and initialisation year
dat_anal [year < ini_yr, maps_pnt := 0]

# calculate incremental cost for RI
dat_anal [, `:=` (inc_cost_del_RI1 = (inc_covMCV1/0.01)*((cov2020_MCV1 < 0.8)*cost_inc_del_RI_lcov +
                                                           (cov2020_MCV1 >= 0.8)*cost_inc_del_RI_hcov) +
                    floor(inc_covMCV1/0.05)*cost_inc_del_RI_per*cost_del_RI,
                  inc_cost_del_RI2 = (inc_covMCV2/0.01)*((cov2020_MCV2 < 0.8)*cost_inc_del_RI_lcov +
                                                           (cov2020_MCV2 >= 0.8)*cost_inc_del_RI_hcov) +
                    floor(inc_covMCV2/0.05)*cost_inc_del_RI_per*cost_del_RI)]

# calculate costs except for MAP procurement costs
dat_anal [, `:=` (all_cost_tx = cases*cost_tx,
                  all_cost_vac_syr = (1-maps_pnt)*((doseRI1+doseRI2)*cost_syr_RI +
                                                     doseSIAc*cost_syr_SIA),
                  all_cost_del_syr = (1-maps_pnt)*(doseRI1*(cost_del_RI+inc_cost_del_RI1) +
                                                     doseRI2*(cost_del_RI+inc_cost_del_RI2) +
                                                     doseSIAc*cost_del_SIA),
                  all_cost_del_maps = maps_pnt*doseRI1*(cost_del_RI+inc_cost_del_RI1) +
                                      doseSIAde1*(cost_del_RI+inc_cost_del_RI1) +
                                      maps_pnt*doseRI2*(cost_del_RI+inc_cost_del_RI2) +
                                      doseSIAde2*(cost_del_RI+inc_cost_del_RI2) +
                                      (maps_pnt*doseSIAc+doseSIAf)*cost_del_SIA,
                  all_dose_maps = maps_pnt*(doseRI1+doseRI2+doseSIAc) + (doseSIAde1+doseSIAde2+doseSIAf))]
dat_anal_full <- copy (dat_anal)

# add income-level data for analysis
dat_anal <- dat_anal_full [, .(country, country_name, income_g, year, comp,
                               dalys, all_dose_maps, all_cost_tx,
                               all_cost_vac_syr, all_cost_del_syr, all_cost_del_maps)]
dat_anal [, costs_exc_procuremaps := all_cost_vac_syr + all_cost_del_syr +  # costs except maps procurement cost
                                     all_cost_tx + all_cost_del_maps]
dat_anal <- dat_anal [, .SD, .SDcols = !c("all_cost_vac_syr", "all_cost_del_syr",
                                          "all_cost_tx", "all_cost_del_maps")]
dat_anal_income <- dat_anal [, lapply(.SD, sum),
                             .SDcols = dalys:costs_exc_procuremaps,
                             by = income_g:comp]
dat_anal <- rbind (dat_anal,
                   dat_anal_income [, `:=` (country = income_g,
                                            country_name = NA)])


# ------------------------------------------------------------------------------
## calculate the maximum prices for introducing MR-MAPs to be cost-effective
# ------------------------------------------------------------------------------
# assume equal discounting rates
discrates <- 1/((1+0.03)^(c(2030:2040)-2020))

# create a function to return ICER given different MR-MAP prices
get_ICER_mapsprice <- function (iso3,
                                scn_eval, # scenario for evaluation
                                scn_base, # baseline scenario
                                price_maps){

  # calculate averted DALYs and incremental costs
  dat_ctry <- dat_anal [country == iso3 & comp %in% c(scn_eval, scn_base)]
  dat_ctry [, costs := costs_exc_procuremaps + all_dose_maps*price_maps]
  avt_ctry <- setDT (pivot_wider (dat_ctry [, .(country, country_name, income_g, year, comp, dalys, costs)],
                                  values_from = c(dalys, costs),
                                  names_from = comp))
  avt_ctry [, `:=` (avt_dalys = get(paste0("dalys_",scn_base)) - get(paste0("dalys_",scn_eval)),
                    inc_costs = get(paste0("costs_",scn_eval)) - get(paste0("costs_",scn_base)))]

  # calculate ICER with 3% equal discounting rates to both health and cost measurements
  setorder (avt_ctry, year) # sort data by year
  total_inc_costs_disc <- sum(avt_ctry$inc_costs*discrates)
  total_avt_dalys_disc <- sum(avt_ctry$avt_dalys*discrates)
  return (list (total_inc_costs_disc = total_inc_costs_disc,
                total_avt_dalys_disc = total_avt_dalys_disc,
                icer = ifelse (total_avt_dalys_disc < 0, NA, total_inc_costs_disc/total_avt_dalys_disc)))
}

# create an empty table to store the results
dat_mapsprice <- dat_anal [year == 2030 & comp %in% scns[c(2,3,5,6)],
                           .(country, country_name, income_g, comp)]
dat_mapsprice [, `:=` (maps_price = as.numeric(NA), obj_diff = as.numeric(NA),
                       icer = as.numeric(NA), icer_threshold = as.numeric(NA))]

# find the threshold prices
for (ictry in unique (dat_anal$country)){
  c_icer_oppcost <- ifelse (ictry %in% anl_countries,
                            dat_oppcost_ctry [country_code == ictry, oppcost],
                            dat_oppcost_income [income_g == ictry, wt_oppcost])
  for (iscn_base in scns[c(1,4)]){
    if (iscn_base == scns[1]){
      scns_to_eval <- scns[2:3]
      } else {
      scns_to_eval <- scns[5:6]
    }
    for (iscn_eval in scns_to_eval){

      # check cost-effectiveness if MR-MAPs procurement is provided free
      icer_0 <- get_ICER_mapsprice (ictry, iscn_eval, iscn_base, 0)$icer
      if (is.na(icer_0) || (icer_0 > c_icer_oppcost)){
        dat_mapsprice [country == ictry & comp == iscn_eval,
                       `:=`(maps_price = NA, obj_diff =  NA, icer = NA,
                            icer_threshold = c_icer_oppcost)]
      } else {

        # find price threshold if (1) positive averted DALYs & (2) icer < threshold when provided free
        obj_fun <- function (price_maps) {
          icer_res <- get_ICER_mapsprice (ictry, iscn_eval, iscn_base, price_maps)$icer
          return (abs(icer_res - c_icer_oppcost))
        }
        optim_res <- optimize (obj_fun, c(1e5,0), tol = 1e-6)
        dat_mapsprice [country == ictry & comp == iscn_eval,
                       `:=`(maps_price = optim_res$minimum,
                            obj_diff =  optim_res$objective,
                            icer = get_ICER_mapsprice (ictry, iscn_eval, iscn_base, optim_res$minimum)$icer,
                            icer_threshold = c_icer_oppcost)]
      }
    }
  }
}


# ------------------------------------------------------------------------------
## Table 5 - MR-MAP price thresholds
# ------------------------------------------------------------------------------
# income-level price thresholds
dat_mapsprice [is.na(country_name)]

# ranges of country-level price thresholds
dat_mapsprice [!is.na(country_name) & !is.na(maps_price),
               .(min_price = min(maps_price), max_price = max(maps_price), med_price = median(maps_price)),
               by = .(comp, income_g)]
table (dat_mapsprice [!is.na(country_name) & is.na(maps_price), .(income_g, comp)])
dat_mapsprice [!is.na(country_name) & is.na(maps_price)]

# check countries have no net health benefits even if MR-MAPs procurement is free
ce_neg_ctry <- unique (dat_mapsprice [is.na(maps_price), country])
total_burden <- annual_burden [, .(total_case = sum(cases)), by = .(country, comp)]
ave_burden <- total_burden [, .(mean_case = mean(total_case)), by = .(country)]
sum (ave_burden [country %in% ce_neg_ctry, mean_case]) / sum (ave_burden [, mean_case]) # 6.3% of cases between 2030-2040
sum (dat_input [comp == "high-base" & country_code %in% ce_neg_ctry, sumcase_17to19]) /
  sum (dat_input [comp == "high-base", sumcase_17to19]) # 18% of 2017-2019 reported cases

