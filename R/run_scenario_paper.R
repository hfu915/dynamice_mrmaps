# run_scenario_paper
# execute dynaMICE model for MR-MAPs scenarios
# update: 2022/11/24

rm(list = ls())

library (data.table)
library (stringr)

# rcpp codes
library (Rcpp)
sourceCpp ("Rcpp/rcpp_spinup_maps.cpp")
sourceCpp ("Rcpp/rcpp_vaccine_oney_maps.cpp")

# load data
load (file = "data/data_pop_maps.rda")
load (file = "data/data_cfr_portnoy_maps.rda")
load (file = "data/data_contact_syn_maps.rda")
load (file = "data/data_r0_maps_paper.rda")     # assume fixed R0 by income group
load (file = "data/data_timeliness_maps.rda")   # assume perfect timeliness
load (file = "data/data_lexp_remain_maps.rda")

load (file = "inputs/dat_income.rds")           # input income data for 102 countries

# prepare functions
source ("R/logs.R")
source ("R/utils_maps.R")
source ("R/functions_rcpp_maps.R")

var <- list (
  vaccine_coverage_folder  = "vac_coverage_maps/",
  burden_estimate_folder   = "burden_estimate/",
  countries                = dat_income [income_g != "High income", country_code]  # 90 countries for evaluation
)

dir.create (file.path (paste0 (getwd(), "/", var$burden_estimate_folder)), recursive = T)

# prepare coverage inputs
vaccine_strategies <- c("high-base", "high-maps", "high-accl",
                         "low-base",  "low-maps",  "low-accl")

# run model
for (index in 1:length(vaccine_strategies)){

  scenario_name  <- vaccine_strategies [index]
  print (scenario_name)
  scenario_number <- sprintf ("scenario%02d", index)

  # run model and estimate cases
  burden_estimate <- runScenario_rcpp (
    vaccine_coverage_folder    = var$vaccine_coverage_folder,
    scenario_name              = scenario_name,
    save_scenario              = scenario_number,
    burden_estimate_folder     = var$burden_estimate_folder,
    countries                  = var$countries,
    vaccination                = 2,  # with both MCV1 and MCV2
    using_sia                  = 1,  # with campaign activities
    contact_mat                = "syn",
    sim_years                  = 1980:2040
  )

  burden_estimate_file <- paste0 ("burden_", scenario_name, ".csv")

  # output the summary csv file
  # estimate deaths and DALYs by Portnoy's CFR methods
  output_burden_estimate (save_scenario          = sprintf ("scenario%02d", index),
                          foldername             = paste0 ("20221124", "_v2_s1_deter"),
                          sim_years              = 1980:2040,
                          cfr_option             = "Portnoy",
                          burden_estimate_file   = burden_estimate_file,
                          burden_estimate_folder = var$burden_estimate_folder,
                          vimc_scenario          = "campaign-default",  # scenario name used used in Portnoy's CFRs
                          portnoy_scenario       = "s6")                # constant CFR after 2018
}

