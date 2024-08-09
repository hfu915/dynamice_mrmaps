# output_results.R
# get tables and figures for the manuscript
# include both health impact and cost-effectiveness analysis
# update: 2024/08/09

library(data.table)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)


rm (list = ls())
# load input information:
# (1) data_input: MAP intro year, unit costs, coverage, reported cases, funding groups, WHO region, Gavi-eligibility
# (2) data_oppcost_income: income-level health opportunity costs
# (3) data_oppcost_ctry: country-level health opportunity costs
# (4) data_pop: UNWPP 2019 population size by age and year
# (5) data_income: country list by income group

load (file = "data/data_input.rds")
load (file = "data/data_oppcost_income.rds")
load (file = "data/data_oppcost_ctry.rds")
load (file = "data/data_pop_maps.rda")
load (file = "data/data_income.rds")

# ------------------------------------------------------------------------------
## exclude countries with small number of cases over 2027-2029
# ------------------------------------------------------------------------------
# load data in the lower coverage, without MR-MAPs scenario
est_low_base <- fread ("burden_estimate/burden_low-base_Portnoy.csv")
cum_low_base <- est_low_base [, .(cases = sum(cases), pops = sum(cohort_size)),
                              by = .(year, country)]
cum_low_base [, incd := cases/pops]
tmptb1  <- table (cum_low_base [year %in% 2027:2029 & incd < 1e-6, country])
names (tmptb1 [tmptb1 >= 3])
# [1]  "ALB" "BGD" "BWA" "CHN" "CPV" "DJI" "ERI" "FJI" "GEO" "GHA"
# [11] "GIN" "GNB" "JOR" "KAZ" "KGZ" "MAR" "MUS" "MWI" "MYS" "NER"
# [21] "ROU" "RUS" "RWA" "SDN" "SOM" "SSD" "THA" "TON" "TUR" "UGA" "VNM"
tmptb2  <- table (cum_low_base [year %in% 2027:2029 & cases < 5, country])
low_burden_ctries <- names (tmptb2 [tmptb2 >= 3]) # all have low incidence rates
low_burden_ctries
# [1]  "ALB" "BWA" "CPV" "DJI" "ERI" "FJI" "GEO" "GIN" "GNB" "JOR"
# [11] "KAZ" "KGZ" "MAR" "MUS" "MYS" "RWA" "SSD" "THA" "TON" "VUT"
tmptb3  <- setDT (pivot_wider (cum_low_base[year %in% 2027:2029,
                                        .(country, year, cases)],
                               names_from = year, values_from = cases))
tmptb3 [, sumcases := `2027` + `2028` + `2029`]
tmptb3 [sumcases <= 10, country]
# [1]  "ALB" "BWA" "CPV" "DJI" "ERI" "FJI" "GEO" "GNB" "JOR" "KAZ"
# [11] "KGZ" "MAR" "MUS" "MYS" "RWA" "SSD" "TON" "VUT"


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
anl_countries <- setdiff (data_income [income_g != "High income", country_code],
                          low_burden_ctries) # 67
table (data_income [country_code %in% anl_countries, income_g])
# Low income          Lower middle income   Upper middle income
# 20                  35                    15

# number of cases between 2030-2040
file_case  <- NULL
for (scname in scns){
  scn_case <- fread (paste0 ("burden_estimate/burden_", scname, "_Portnoy.csv"))
  scn_case [, comp := scname]
  file_case <- rbind (file_case,
                      scn_case [country %in% anl_countries & year %in% 2030:2040])
}
annual_burden <- file_case [, lapply(.SD, sum),
                            .SDcols = cases:doseSIAf2,
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
dat_anal <- annual_burden [data_input [country_code %in% anl_countries],
                           on = .(country = country_code, comp = comp)]
setnames (x = dat_anal, old = "i.country", new = "country_name")

# adjust market penetration rate based on scenario and introduction year
dat_anal [year < ini_yr, maps_pnt := 0]

# calculate incremental cost for RI
dat_anal [, `:=` (inc_cost_del_RI1 = (inc_covMCV1/0.01)*((cov2020_MCV1 < 0.8)*cost_inc_del_RI_lcov +
                                                           (cov2020_MCV1 >= 0.8)*cost_inc_del_RI_hcov) +
                    floor(inc_covMCV1/0.05)*cost_inc_del_RI_per*cost_del_RI,
                  inc_cost_del_RI2 = (inc_covMCV2/0.01)*((cov2020_MCV2 < 0.8)*cost_inc_del_RI_lcov +
                                                           (cov2020_MCV2 >= 0.8)*cost_inc_del_RI_hcov) +
                    floor(inc_covMCV2/0.05)*cost_inc_del_RI_per*cost_del_RI)]
dim (dat_anal) # 4620 (6 scenarios, 70 countries, 11 years)

# sum total doses and costs
# market penetration rate applied to strategies A, B, and C
# costs are decomposed for further analysis
# generate costs by delivery costs
get_cost_del <- function (
    rel_del_maps = 1,      # relative reduction of delivery cost for MR-MAPs compared to N&Ss
    diff_del_maps_RI = 0,  # difference in delivery cost for MR-MAPs compared to N&Ss (RI)
    diff_del_maps_SIA = 0  # difference in delivery cost for MR-MAPs compared to N&Ss (SIA)
){
  dat_anal [, `:=` (total_doses = doseRI1+doseRI2 + doseSIAc1+doseSIAc2 +
                      doseSIAde1+doseSIAde2 + doseSIAf1+doseSIAf2,
                    all_cost_tx = cases*cost_tx,
                    all_cost_vac_syr = (1-maps_pnt)*((doseRI1+doseRI2)*cost_syr_RI +
                                                       (doseSIAc1+doseSIAc2)*cost_syr_SIA),
                    all_cost_del_syr = (1-maps_pnt)*(doseRI1*(cost_del_RI+inc_cost_del_RI1) +
                                                       doseRI2*(cost_del_RI+inc_cost_del_RI2) +
                                                       (doseSIAc1+doseSIAc2)*cost_del_SIA),
                    all_cost_vac_maps_lb = maps_pnt*(doseRI1+doseRI2+doseSIAc1+doseSIAc2)*cost_map_lb +
                      (doseSIAde1+doseSIAde2+doseSIAf1+doseSIAf2)*cost_map_lb,
                    all_cost_vac_maps_ub = maps_pnt*(doseRI1+doseRI2+doseSIAc1+doseSIAc2)*cost_map_ub +
                      (doseSIAde1+doseSIAde2+doseSIAf1+doseSIAf2)*cost_map_ub,
                    all_cost_del_maps = ((maps_pnt*doseRI1*(cost_del_RI+inc_cost_del_RI1+diff_del_maps_RI) +
                                            doseSIAde1*(cost_del_RI+inc_cost_del_RI1+diff_del_maps_RI) +
                                            maps_pnt*doseRI2*(cost_del_RI+inc_cost_del_RI2+diff_del_maps_RI) +
                                            doseSIAde2*(cost_del_RI+inc_cost_del_RI2+diff_del_maps_RI) +
                                            maps_pnt*(doseSIAc1+doseSIAc2)*(cost_del_SIA+diff_del_maps_SIA) +
                                            (doseSIAf1+doseSIAf2)*(cost_del_SIA+diff_del_maps_SIA))*rel_del_maps))]
  dat_anal [, `:=` (total_costs_lb = all_cost_tx +
                      all_cost_vac_syr + all_cost_del_syr +
                      all_cost_vac_maps_lb + all_cost_del_maps,
                    total_costs_ub = all_cost_tx +
                      all_cost_vac_syr + all_cost_del_syr +
                      all_cost_vac_maps_ub + all_cost_del_maps)]

  # get results by income level and global level
  dat_country <- copy (dat_anal [, .(comp, income_g, country_name, country, year,
                                     cohort_size, cases, deaths, dalys,
                                     total_doses, all_cost_tx,
                                     all_cost_vac_syr, all_cost_del_syr,
                                     all_cost_vac_maps_lb, all_cost_vac_maps_ub, all_cost_del_maps,
                                     total_costs_lb, total_costs_ub)]) [, lapply (.SD, sum),
                                                                        .SDcols = cohort_size:total_costs_ub,
                                                                        by = c("comp", "income_g", "country_name", "country", "year")]
  dat_income <- copy (dat_anal [, .(comp, income_g, country_name, country, year,
                                    cohort_size, cases, deaths, dalys,
                                    total_doses, all_cost_tx,
                                    all_cost_vac_syr, all_cost_del_syr,
                                    all_cost_vac_maps_lb, all_cost_vac_maps_ub, all_cost_del_maps,
                                    total_costs_lb, total_costs_ub)]) [, lapply (.SD, sum),
                                                                       .SDcols = cohort_size:total_costs_ub,
                                                                       by = c("comp", "income_g", "year")]
  dat_income [, ':=' (country = income_g, country_name = income_g)]
  dat_global <- copy (dat_anal [, .(comp, income_g, country_name, country, year,
                                    cohort_size, cases, deaths, dalys,
                                    total_doses, all_cost_tx,
                                    all_cost_vac_syr, all_cost_del_syr,
                                    all_cost_vac_maps_lb, all_cost_vac_maps_ub, all_cost_del_maps,
                                    total_costs_lb, total_costs_ub)]) [, lapply (.SD, sum),
                                                                       .SDcols = cohort_size:total_costs_ub,
                                                                       by = c("comp", "year")]
  dat_global [, ':=' (country = "Global", country_name = "Global", income_g = "Global")]

  dat_all <- rbind (dat_country,
                    setcolorder (dat_income, names(dat_country)),
                    setcolorder (dat_global, names(dat_country)))
  return(dat_all)
}

dat_all <- get_cost_del(1,0,0)


# ------------------------------------------------------------------------------
## calculate cost-effectiveness measures
# ------------------------------------------------------------------------------
# assign discounting rates
get_avt_cea <- function (caldat, base_scn, dr_cost, dr_eff, dr_refyr){
  dat_anal_base <- caldat [comp == base_scn]
  dat_anal_intr <- caldat [dat_anal_base[, c("year", "country", "cases", "deaths", "dalys",
                                             "total_doses", "total_costs_lb", "total_costs_ub")],
                           on = .(country = country, year = year)]
  setnames (x = dat_anal_intr,
            old = c("comp", "i.cases", "i.deaths", "i.dalys",
                    "i.total_doses", "i.total_costs_lb", "i.total_costs_ub"),
            new = c("scenario", "cases_base", "deaths_base", "dalys_base",
                    "total_doses_base", "total_costs_base_lb", "total_costs_base_ub"))

  dat_disc <- data.table (year = 2030:2040,
                          disc_cost = 1/((1+dr_cost)^(c(2030:2040)-dr_refyr)),
                          disc_eff = 1/((1+dr_eff)^(c(2030:2040)-dr_refyr)))
  dat_anal_intr  <- dat_anal_intr [dat_disc, on = .(year = year)]

  dat_anal_intr [, `:=` (avt_dalys = (dalys_base - dalys)*disc_eff,
                         avt_cases = cases_base - cases,
                         avt_deaths = deaths_base - deaths,
                         inc_doses = total_doses - total_doses_base,
                         inc_cost_lb = (total_costs_lb - total_costs_base_lb)*disc_cost,
                         inc_cost_ub = (total_costs_ub - total_costs_base_ub)*disc_cost)]

  dat_anal_intr_sum <-  dat_anal_intr [, lapply(.SD, sum),
                                       .SDcols = !c("year", "disc_cost", "disc_eff"),
                                       by = c("country", "country_name", "income_g", "scenario")]

  # calculate (1) proportion of case reduction
  #           (2) incremental cost-effectiveness ratio
  #           (3) number of doses needed to be given
  dat_anal_intr_sum [, `:=` (prred_cases   = avt_cases/cases_base,
                             prred_deaths  = avt_deaths/deaths_base,
                             prred_dalys   = avt_dalys/dalys_base,
                             nnd_case      = inc_doses/avt_cases,
                             icer_lb       = inc_cost_lb/avt_dalys,
                             icer_ub       = inc_cost_ub/avt_dalys,
                             icer_case_lb  = inc_cost_lb/avt_cases,
                             icer_case_ub  = inc_cost_ub/avt_cases,
                             icer_death_lb = inc_cost_lb/avt_deaths,
                             icer_death_ub = inc_cost_ub/avt_deaths)]

  return (dat_anal_intr_sum)
}

dat_cea_high <- get_avt_cea (dat_all [comp %in% scns[1:3]], scns[1], 0, 0, 2020)
dat_cea_low  <- get_avt_cea (dat_all [comp %in% scns[4:6]], scns[4], 0, 0, 2020)


# ------------------------------------------------------------------------------
# set up plotting elements
# ------------------------------------------------------------------------------
plot_scn_names <- c("Without MR-MAPs\n (higher coverage)",
                    "Sequential intro\n (higher coverage)",
                    "Accelerated intro\n (higher coverage)",
                    "Without MR-MAPs\n (lower coverage)",
                    "Sequential intro\n (lower coverage)",
                    "Accelerated intro\n (lower coverage)")
strategy_names <- c("A: Routine MCV1, Surviving children, 9m",
                    "B: Routine MCV2, Surviving children, 16.5m",
                    "C: SIA, Varying age groups",  # after 2030, all are 9-59m
                    "D: Routine MCV1, HTR & MOV, 1-2y",
                    "E: Routine MCV2, HTR & MOV, 1-2y",
                    "F: One-time catch-up, HTR & MOV, 2-15y")
income_names <- c("Low income", "Lower middle income", "Upper middle income")
plot_colours <- c("#1f78b4", "#33a02c", "#e31a1c",
                  "#a6cee3", "#b2df8a", "#fb9a99")
# blank plot
blankPlot <- ggplot() + geom_blank(aes(1,1)) +
  theme (plot.background = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.text.x = element_blank(),
         axis.text.y = element_blank(),
         axis.ticks = element_blank(),
         axis.line = element_blank())


# ------------------------------------------------------------------------------
# Fig 1: coverage for delivery strategies in a selected country
# ------------------------------------------------------------------------------
sel_cols_f1 <- c("strategy", "country_code",  "country", "year", "coverage")
plot_cov_ctry <- function (sel_ctry, show_ctry_name = TRUE){
  file_cov_ctry <- NULL
  for (scname in scns){
    input_scn_cov <-  rbind (fread (paste0 ("vac_coverage_maps/routine_",
                                            scname, ".csv")) [country_code == sel_ctry & year >= 2020, ..sel_cols_f1],
                             fread (paste0 ("vac_coverage_maps/sia_",
                                            scname, ".csv")) [country_code == sel_ctry & year >= 2020, ..sel_cols_f1])
    file_cov_ctry <- rbind (file_cov_ctry, input_scn_cov [, comp := scname])
  }
  file_cov_ctry <- file_cov_ctry [data_input [country_code == sel_ctry, c("comp", "ini_yr")],
                                on = .(comp)]

  # add an empty row for strategy C to generate the panel
  if (length(file_cov_ctry [year %in% 2030:2040 & str_sub(strategy,1,1) == "C", coverage]) == 0){
    file_cov_ctry <- rbind (file_cov_ctry,
                            copy(file_cov_ctry[year == 2040 &
                                                 str_sub(strategy,1,1) == "A"])[, `:=`(strategy = strategy_names[3],
                                                                                       coverage = NA)])
  }
  file_cov_ctry [, year := as.numeric(year)]

  file_cov_ctry [comp %in% scns[1:3] & str_sub(strategy,1,1) %in% c("D","E","F"), year := year - 0.15]
  file_cov_ctry [comp %in% scns[4:6] & str_sub(strategy,1,1) %in% c("D","E","F"), year := year + 0.15]
  file_cov_ctry [, `:=` (scenario = factor (comp, levels = scns, labels = plot_scn_names),
                         strategy = factor (str_sub(strategy,1,1), levels = c("A","B","C","D","E","F"),
                                            labels = strategy_names))]

  # use coverage among target population (instead of general population)
  file_cov_ctry [, adj_cov := coverage]
  file_cov_ctry [str_sub(strategy,1,1) %in% c("D","E"), adj_cov := 0.2]
  file_cov_ctry [str_sub(strategy,1,1) == "F", adj_cov := 0.1]

  plot_f1_cov <- ggplot (data = file_cov_ctry) +
    geom_vline (aes(xintercept = ini_yr, colour = scenario),
                linewidth = 0.75, linetype = 2, show.legend = F) +
    geom_point (aes (x = year, y = adj_cov, colour = scenario, size = scenario)) +
    ylim (0,1) + xlim (2019.5,2040.5) +
    facet_wrap (vars(strategy), ncol = 2, dir = "v") +
    scale_colour_manual ("Scenario", values = plot_colours) +
    scale_size_manual (" ", values = c(rep(c(2.9,1.9,0.8),2))) +
    labs (x = "Year", y = "Coverage",
          title = ifelse (show_ctry_name, file_cov_ctry$country[1], " ")) +
    guides (colour = guide_legend (nrow = 3),
            size = guide_none()) +
    theme_bw () +
    theme (legend.position = "bottom",
           legend.text = element_text(size = 9),
           legend.key.size = unit(0.85, 'cm'))
  return (plot_f1_cov)
}

ggsave ("outputs/paper_fig_cov-COD.pdf", plot_cov_ctry ("COD", FALSE),
        height = 16, width = 17, units = "cm" )

pdf ("outputs/paper_fig_cov-all.pdf", height = 6, width = 7)
for (plt_ctry in anl_countries) {
  print (plot_cov_ctry (plt_ctry, TRUE))
}
dev.off()


# ------------------------------------------------------------------------------
## Table 3 -  absolute and relative averted burden
# ------------------------------------------------------------------------------
dat_cea_high [country %in% c(income_names, "Global"),
              .(country, scenario, avt_cases, avt_deaths, avt_dalys,
                prred_cases, prred_deaths, prred_dalys)]
dat_cea_low  [country %in% c(income_names, "Global"),
              .(country, scenario, avt_cases, avt_deaths, avt_dalys,
                prred_cases, prred_deaths, prred_dalys)]

dat_cea_high [avt_cases < 0] # MLI/Sequential & Accelerated
dat_cea_low  [avt_cases < 0] # AGO/Sequential

dat_cea_high [avt_dalys < 0] # MLI/Accelerated
# Note that MLI/Accelerated have negative cases but positive DALYs averted
# a large outbreak would be delayed to occur and mean age of infection would increase


# ------------------------------------------------------------------------------
## Fig 2 - cumulative burden by income
# ------------------------------------------------------------------------------
pltdata_cumburden <- dat_all [country_name %in% income_names, comp:dalys]
pltdata_cumburden <- pltdata_cumburden [, lapply (.SD, sum),
                                        .SDcols = cases:dalys,
                                        by = c("comp", "income_g")]
# burden size reported in the manuscript text
pltdata_cumburden [, .(sum_cases = sum(cases),
                       sum_deaths = sum(deaths),
                       sum_dalys = sum(dalys)), by = "comp"]

pltdata_cumburden <- setDT (pivot_longer (pltdata_cumburden, cols = cases:dalys,
                                          names_to = "measure", values_to = "value"))
pltdata_cumburden [, `:=` (comp = factor (comp, levels = scns, labels = plot_scn_names),
                           measure = factor (measure, levels = c("cases", "deaths", "dalys"),
                                             labels = c("Cases", "Deaths", "DALYs")))]
plt_cumburden <- ggplot(data = pltdata_cumburden,
                        aes(x = comp, y = value/1e6, fill = comp)) +
  geom_col (position = "dodge") +
  facet_grid (cols = vars(income_g), rows = vars(measure), scales = "free") +
  scale_fill_manual ("Scenarios", values = plot_colours) +
  theme_bw () +
  labs (x = " ", y = "Burden estimates (millions)") +
  guides (fill = guide_legend (byrow = TRUE)) +
  theme (legend.text = element_text (size = 10),
         legend.spacing.y = unit (0.25, "cm"),
         legend.key.size = unit (0.7, "cm"),
         strip.text.x = element_text (size = 11),
         strip.text.y = element_text (size = 11),
         plot.margin = unit (c(0.5, 0.5, 0.5, 0.5), "cm"),
         axis.text.x = element_blank (),
         axis.ticks.x = element_blank (),
         axis.title.y = element_text (size = 13, margin = margin (r = 10)))

ggsave("outputs/paper_fig_cumburden-overview.pdf", plt_cumburden,
       height = 6, width = 9)



# ------------------------------------------------------------------------------
## Table C -  total incremental cost
# ------------------------------------------------------------------------------
dat_cea_high [country %in% c(income_names, "Global"),
              .(country, scenario, inc_doses, inc_cost_lb, inc_cost_ub)]
dat_cea_low  [country %in% c(income_names, "Global"),
              .(country, scenario, inc_doses, inc_cost_lb, inc_cost_ub)]


# ------------------------------------------------------------------------------
## Fig 3 - incremental costs breakdown
# ------------------------------------------------------------------------------
price_suffix <- paste0 (" - ", c("Lower","Upper"), " MR-MAP price")
price_names <- list ("Higher" = c(paste0 (scns[2], price_suffix), paste0 (scns[3], price_suffix)),
                     "Lower" = c(paste0 (scns[5], price_suffix), paste0 (scns[6], price_suffix)))

plt_decomp_inccost <- function (in_data, sel_covassum, sel_title, sel_lgdpos){
  plt_data <- in_data [, lapply (.SD, sum),
                       .SDcols = cohort_size:total_costs_ub,
                       by = c("country", "country_name", "income_g", "comp")]

  scn_base_num <- ifelse (sel_covassum == "Higher", 1, 4)
  plt_data_base <- plt_data [comp == scns [scn_base_num], !(cohort_size:total_doses)]
  plt_data <- plt_data [comp %in% scns [scn_base_num+c(1:2)]]
  plt_data <- plt_data [plt_data_base, on = .(country, country_name, income_g)]
  plt_data [, `:=` (inc_cost_tx          = all_cost_tx - i.all_cost_tx,
                    inc_cost_vac_syr     = all_cost_vac_syr - i.all_cost_vac_syr,
                    inc_cost_del_syr     = all_cost_del_syr - i.all_cost_del_syr,
                    inc_cost_vac_maps_lb = all_cost_vac_maps_lb - i.all_cost_vac_maps_lb,
                    inc_cost_vac_maps_ub = all_cost_vac_maps_ub - i.all_cost_vac_maps_ub,
                    inc_cost_del_maps    = all_cost_del_maps - i.all_cost_del_maps
  )]

  plt_data <- setDT (pivot_longer (plt_data [, !(cohort_size:i.total_costs_ub)],
                                   cols = inc_cost_tx:inc_cost_del_maps,
                                   names_to = "cost_type",
                                   values_to = "inc_cost"))

  dat_plot_decost <- plt_data [country %in% income_names]
  dat_plot_decost <- rbind (copy (dat_plot_decost [cost_type != "inc_cost_vac_maps_ub"]) [, comp := paste0 (comp, " - Lower MR-MAP price")],
                            copy (dat_plot_decost [cost_type != "inc_cost_vac_maps_lb"]) [, comp := paste0 (comp, " - Upper MR-MAP price")])

  dat_plot_decost [, `:=` (comp = factor (comp, levels = price_names[[sel_covassum]],
                                          labels = c("Sequential intro /\n Lower price",
                                                     "Sequential intro /\n Upper price",
                                                     "Accelerated intro /\n Lower price",
                                                     "Accelerated intro /\n Upper price")),
                           cost_type = factor (cost_type,
                                               levels = c("inc_cost_tx",
                                                          "inc_cost_vac_syr", "inc_cost_del_syr",
                                                          "inc_cost_vac_maps_lb", "inc_cost_vac_maps_ub",
                                                          "inc_cost_del_maps"),
                                               labels = c("Treatment of measles illness",
                                                          "Vaccine procurement (N&S)", "Vaccine delivery (N&S)",
                                                          "Vaccine procurement (MR-MAP, lower price)", "Vaccine procurement (MR-MAP, upper price)",
                                                          "Vaccine delivery (MR-MAP)")),
                           income_g = factor (income_g, levels = income_names))]

  dat_plot_decost_total <- dat_plot_decost [, .(inc_cost_total = sum(inc_cost)),
                                            by = country:comp]

  ggplot (data = dat_plot_decost,
          aes (x = comp, y = inc_cost/1e6, group = comp)) +
    geom_col (position = "stack", width = 0.8, aes (fill = cost_type)) +
    geom_hline (yintercept = 0, color = "darkgrey") +
    geom_point (data = dat_plot_decost_total,
                aes(y = inc_cost_total/1e6, group = comp),
                shape = 21, stroke = 2.5, size = 2.5) +
    coord_flip() +
    facet_wrap (vars(income_g), nrow = 1, scales = "free_x") +
    scale_fill_manual ("Cost type", values = c("#fdb863", "#7b3294", "#c2a5cf",
                                               "#31a354", "#006d2c", "#bae4b3")) +
    labs (y = "Incremental cost (millions)", x = " ",
          title = paste0 (sel_covassum, " coverage projection", sel_title)) +
    theme_bw () +
    theme (legend.position = sel_lgdpos,
           legend.direction = "horizontal",
           plot.margin = unit (c(0.5, 0.5, 0.25, 0.25), "cm"),
           strip.text.x = element_text (size = 10),
           axis.text.x = element_text (size = 10.5),
           legend.text = element_text (size = 9.5)) +
    guides (fill = guide_legend (nrow = 2, byrow = T))
}

ggsave ("outputs/paper_fig_inccost-breakdown.pdf",
        ggarrange (plt_decomp_inccost (dat_all, "Higher", "", "none"),
                   plt_decomp_inccost (dat_all, "Lower", "", c(0.45, -0.6)),
                   blankPlot, heights = c(3,3,1), nrow = 3,
                   labels = c("A", "B", " "), label.x = 0.02, label.y = 0.95),
        height = 6.5, width = 10)


# ------------------------------------------------------------------------------
## Table 4 & Table D - number of countries that consider MR-MAPs cost-effective
# ------------------------------------------------------------------------------
# generate ICERs for different discounting approaches
dat_cea_high_dcost <- get_avt_cea (dat_all [comp %in% scns[1:3]], scns[1], 0.03, 0, 2020)
dat_cea_low_dcost  <- get_avt_cea (dat_all [comp %in% scns[4:6]], scns[4], 0.03, 0, 2020)
dat_cea_high_dcosteff <- get_avt_cea (dat_all [comp %in% scns[1:3]], scns[1], 0.03, 0.03, 2020)
dat_cea_low_dcosteff  <- get_avt_cea (dat_all [comp %in% scns[4:6]], scns[4], 0.03, 0.03, 2020)

dat_cea_high_dcosteff [avt_dalys < 0] # NULL
dat_cea_low_dcosteff  [avt_dalys < 0] # AGO/Sequential
# Note that MLI/Sequential no longer have negative DALYs averted after discounting

dat_cea_all <- rbind (copy (dat_cea_high_dcost [scenario != scns[1]])[, `:=` (covassum = "Higher", discount = "differential")],
                      copy (dat_cea_low_dcost [scenario != scns[4]])[, `:=` (covassum = "Lower", discount = "differential")],
                      copy (dat_cea_high_dcosteff [scenario != scns[1]])[, `:=` (covassum = "Higher", discount = "equal")],
                      copy (dat_cea_low_dcosteff [scenario != scns[4]])[, `:=` (covassum = "Lower", discount = "equal")])

# country-level thresholds
dat_cea_all_ctry <- copy (data_oppcost_ctry [, c("country_code", "oppcost")]) [
  dat_cea_all [, c("covassum", "discount", "country", "country_name", "income_g",
                   "scenario", "avt_dalys", "inc_cost_lb", "inc_cost_ub",
                   "icer_lb", "icer_ub")],
  on = .(country_code = country)]
dat_cea_all_cost <- setDT (pivot_longer (dat_cea_all_ctry [, !c("icer_lb", "icer_ub", "avt_dalys",
                                                                "income_g", "country_name", "oppcost")],
                                         cols = inc_cost_lb:inc_cost_ub,
                                         values_to = "inc_cost",
                                         names_to = "map_price",
                                         names_pattern = "inc_cost_(.*)"))
dat_cea_all_ctry <- setDT (pivot_longer (dat_cea_all_ctry [, !c("inc_cost_lb", "inc_cost_ub")],
                                         cols = icer_lb:icer_ub,
                                         values_to = "icer",
                                         names_to = "map_price",
                                         names_pattern = "icer_(.*)"))
dat_cea_all_ctry <- dat_cea_all_ctry [dat_cea_all_cost,
                                      on = c("covassum", "discount", "country_code", "scenario", "map_price")]
dat_cea_all_ctry [,`:=` (yes_ce = ifelse (icer <= oppcost, 1, 0),
                         income_g = factor (income_g, levels = income_names),
                         scenario = factor (scenario, levels = scns[c(2,3,5,6)], labels = plot_scn_names[c(2,3,5,6)]))]

# numbers of countries where introducing MR-MAPs would be cost-effectiveness
table (dat_cea_all_ctry [discount == "differential" & map_price == "lb" & yes_ce == 1,
                         .(income_g, scenario, map_price)])
table (dat_cea_all_ctry [discount == "differential" & map_price == "ub" & yes_ce == 1,
                         .(income_g, scenario, map_price)])
table (dat_cea_all_ctry [discount == "equal" & map_price == "lb" & yes_ce == 1,
                         .(income_g, scenario, map_price)])
table (dat_cea_all_ctry [discount == "equal" & map_price == "ub" & yes_ce == 1,
                         .(income_g, scenario, map_price)])

# income-level estimates
tab_icer_cols <- c("country", "scenario", "covassum", "discount",
                   "inc_cost_lb", "inc_cost_ub", "icer_lb", "icer_ub")
dat_cea_all_income <- dat_cea_all [country %in% c(income_names, "Global"), ..tab_icer_cols]
dat_cea_all_income <- setDT (pivot_longer (dat_cea_all_income [, !c("inc_cost_lb", "inc_cost_ub")],
                                           cols = icer_lb:icer_ub,
                                           values_to = "icer",
                                           names_to = "map_price",
                                           names_pattern = "icer_(.*)"))
dat_cea_all_income <- dat_cea_all_income [data_oppcost_income[, c("income_g", "wt_oppcost")],
                                          on = .(country = income_g)]
dat_cea_all_income [, yes_ce := ifelse (icer <= wt_oppcost, 1, 0)]


# ------------------------------------------------------------------------------
## Fig 4 - income-level ICERs
# ------------------------------------------------------------------------------
pltdata_icers_income <- copy (dat_cea_all_income)

# transformed to 'special' log scale
pltdata_icers_income [, `:=` (map_price = factor (map_price, levels = c("lb", "ub"),
                                                  labels = c("Lower price", "Upper price")),
                              scenario = factor (scenario, levels = scns[c(2,3,5,6)],
                                                 labels = c("Higher coverage /\n Sequential intro",
                                                            "Higher coverage /\n Accelerated intro",
                                                            "Lower coverage /\n Sequential intro",
                                                            "Lower coverage /\n Accelerated intro")),
                              discount = factor (discount, levels = c("equal", "differential"),
                                                 labels = c("Equal discounting", "Differential discounting")))]
setorder (pltdata_icers_income, discount, country, scenario, map_price)

pltdata_icers_income [, .(icer_lb = min(icer), icer_ub = max(icer)), by = "country"]
 #             country    icer_lb  icer_ub
 #          Low income   10.57211 1846.441
 # Lower middle income -108.09814 1100.733
 # Upper middle income -133.66487 1205.791

pltdata_icers_income [, .(icer_lb = min(icer), icer_ub = max(icer)), by = "map_price"]
#   map_price    icer_lb   icer_ub
# Lower price -133.66487  625.7684
# Upper price   10.06505 1846.4414

plot_icer_income <- function (sel_discount, lgd_pos){
  ggplot (pltdata_icers_income [discount == sel_discount],
          aes(x = scenario, y = icer, colour = scenario)) +
    facet_wrap (vars(country)) +#, scales = "free_y"
    scale_colour_manual ("Scenario", values = plot_colours[c(2,3,5,6)], guide = "none") +
    geom_jitter (aes (shape = map_price), size = 2, stroke = 0.9, width = 0.1) +
    scale_shape_manual ("MR-MAP price", values = c(1, 17)) +
    geom_hline (aes(yintercept = wt_oppcost), colour = "grey60", linetype = 2, linewidth = 0.8) +
      labs (x = "", y = "Cost per DALY averted (USD)", title = sel_discount) +
    theme_bw () +
    theme (legend.position = lgd_pos,
           legend.direction = "horizontal",
           panel.grid.minor = element_blank(),
           plot.margin = unit (c(0.2, 0.25, 0.35, 0.4), "cm"),
           strip.text.x = element_text (size = 12),
           strip.text.y = element_text (size = 12),
           axis.text.y = element_text (size = 10),
           axis.text.x = element_text (size = 10, angle = 60, vjust = 0.55),
           legend.text = element_text (size = 10.5),
           legend.title = element_text (size = 10.5))
}
ggsave ("outputs/paper_fig_icer_income.pdf",
        ggarrange (plot_icer_income ("Equal discounting", "none"),
                   plot_icer_income ("Differential discounting", c(0.5, -0.8)),
                   blankPlot, heights = c(3,3,0.3), nrow = 3,
                   labels = c("A", "B", " "), label.x = 0.01, label.y = 0.993),
        height = 8, width = 10)


# ------------------------------------------------------------------------------
## Fig B - Country-specific net health benefits (equal discounting)
# ------------------------------------------------------------------------------
pltdat_nhb <- copy (dat_cea_all_ctry [discount == "equal"])
pltdat_nhb [, intro := factor (scenario, levels = plot_scn_names[c(2,3,5,6)],
                               labels = c("maps", "accl", "maps", "accl"))]
pltdat_nhb [, `:=` (scenario2 = factor (paste0 (intro, "-", map_price),
                                        levels = c("maps-lb", "accl-lb", "maps-ub", "accl-ub"),
                                        labels = c("Sequential intro/\nLower MR-MAP price",
                                                   "Accelerated intro/\nLower MR-MAP price",
                                                   "Sequential intro/\nUpper MR-MAP price",
                                                   "Accelerated intro/\nUpper MR-MAP price")),
                    map_price = factor (map_price, levels = c("lb", "ub"),
                                        labels = c("Lower price", "Upper price")),
                    nhb = avt_dalys - inc_cost/oppcost)]

pltdat_nhb [, .(med_nhb = median(nhb)),
            by = c("covassum", "income_g", "scenario", "map_price")]

plot_country_nhb <- function (sel_covassum, sel_lgdpos, sel_ytitle){
  ggplot(pltdat_nhb [covassum == sel_covassum],
         aes(x = scenario2, y = nhb/1e6, colour = intro)) +
    facet_wrap (vars(income_g), scales = "free") +
    geom_hline (yintercept = 0, colour = "grey60", linetype = 2, linewidth = 0.8) +
    geom_jitter (aes (shape = map_price), size = 1.2,
                 stroke = 0.9, width = 0.2, alpha = 0.8) +
    scale_colour_manual ("Introduction strategy", values = plot_colours[c(2,3)],
                         labels = c("Sequential", "Accelerated")) +
    scale_shape_manual ("MR-MAP price", values = c(1, 17)) +
    labs (title = paste0 (sel_covassum, " coverage projection"), x = "", y = sel_ytitle) +
    theme_bw () +
    theme (legend.position = sel_lgdpos,
           legend.direction = "horizontal",
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin = unit (c(0.25, 0.25, 0.1, 0.25), "cm"),
           strip.text.x = element_text (size = 10),
           axis.text.x = element_blank(),
           legend.text = element_text (size = 10),
           legend.background = element_blank())
}
ggsave ("outputs/paper_fig_nhb.pdf",
        ggarrange (plot_country_nhb ("Higher", "none", "Net health benefits (millions)"),
                   plot_country_nhb ("Lower", c(0.5, -0.3), "Net health benefits (millions)"),
                   blankPlot,
                   heights = c(3,3,1), ncol = 1,
                   labels = c("A", "B", " "), label.x = 0.01, label.y = 0.99),
        height = 7, width = 10)

