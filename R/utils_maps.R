# utils_maps.R
# Functions for supporting utilities in the DynaMICE model (Dynamic Measles Immunisation Calculation Engine)


# ------------------------------------------------------------------------------
#' Expand the matrix to a different dimension
#
#' This function returns an expanded contact matrix with the specified age
#' structure for model inputs.
#'
#' @param A A matrix to be expanded.
#' @param expand_rows A number of times to repeat each row.
#' @param expand_cols A number of times to repeat each column.
#' @param rescale_rows A logical variable to control whether to rescale the
#' expanded rows.
#' @param rescale_cols A logical variable to control whether to rescale the
#' expanded columns.
#' @return An expanded matrix with rescaling if applicable.
#' @examples
#' expandMatrix (matrix(1:9,3,3), 2, 1, FALSE, FALSE)
expandMatrix <- function (A,
                          expand_rows  = 1,
                          expand_cols  = 1,
                          rescale_rows = F,
                          rescale_cols = F) {

  if(!is.matrix(A)){
    stop("A is not a matrix")
  }

  matvals <- numeric(0)
  rows <- nrow(A)
  cols <- ncol(A)

  for(c in 1:cols) {
    matvals <- c(
      matvals,
      rep(
        A[,c],
        expand_cols,
        each = expand_rows
      )
    )
  }

  B <- matrix (matvals,
               nrow = rows * expand_rows,
               ncol = cols * expand_cols)

  if(rescale_rows & rescale_cols){
    B <- B/(expand_rows*expand_cols)
  } else if(rescale_rows){
    B <- B/expand_rows
  } else if(rescale_cols){
    B <- B/expand_cols
  }

  return (B)

} # end of function -- expandMatrix
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
#' Tailor the data structure for life expectancy by age and year
#'
#' Tailor the \code{\link{data_lexp_remain}} data to the format for processing burden
#' estimates. Linear interpolation between calender years was applied.
#'
#' @param sel_countries ISO-3 codes for countries included for evaluation. If
#' "all", all countries in the original data are selected.
#'
#' @import data.table
#'
#' @examples
#' lexp_remain <- tailor_data_lexp_remain (sel_countries = c("AGO","BGD"))
tailor_data_lexp_remain <- function (sel_countries = "all",
                                     sel_years = 1980:2100){

  # select years included in the dataset for intrapolation
  dat_years <- sel_years [sel_years%%5 == 0]
  if ( min(dat_years) > min(sel_years) ) {dat_years <- c(min(dat_years)-5, dat_years)}
  if ( max(dat_years) < max(sel_years) ) {dat_years <- c(dat_years, max(dat_years)+5)}
  dat_years <- dat_years [dat_years >= 1950 & dat_years <= 2095]
  lexp_remain <- copy (setDT (data_lexp_remain)) [year %in% dat_years]

  if (sel_countries[[1]] != "all") {
    lexp_remain <- lexp_remain [country_code %in% sel_countries]
  }

  # calculate the difference between years
  setorder (lexp_remain, country_code, age_from, year)
  lexp_remain [ , diffy := value - shift(value) , by = .(country_code, age_from)]

  # interpolate a linear trend for between-years
  lexp_remain_yr <- copy (lexp_remain) [year == min(dat_years)]
  for (btwyr in 0:4) {
    dt <- copy (lexp_remain) [year != min(dat_years)]
    dt [, `:=` (year = year - btwyr,
                value = value - btwyr*(diffy/5))]

    lexp_remain_yr <- rbind (lexp_remain_yr, dt)
  }
  lexp_remain_yr <- copy (lexp_remain_yr [year %in% sel_years]) [, diffy := NULL]

  # interpolate a linear trend for between-ages
  setorder (lexp_remain_yr, country_code, year, age_from)
  lexp_remain_yr [ , age_mean := (age_from + age_to)/2]

  lexp_remain_yr_age <- lexp_remain_yr [, .(value = approx(age_mean, value, xout = 0:100)$y),
                                      by = .(country_code, country, year)]
  lexp_remain_yr_age [, age := rep(0:100, length(unique(year))*length(unique(country_code)))]

  # copy values for years before 1950 years after 2095
  lexp_remain_full <- rbind (lexp_remain_yr_age,
                             rbindlist (lapply (min(sel_years):1949, function(iyr) copy (lexp_remain_yr_age [year == 1950])[, year := iyr])),
                             rbindlist (lapply (2096:max(sel_years), function(iyr) copy (lexp_remain_yr_age [year == 2095])[, year := iyr])))

  setorder (lexp_remain_full, country_code, year, age)
  return (lexp_remain_full)

} # end of function -- tailor_data_lexp_remain
# ------------------------------------------------------------------------------

