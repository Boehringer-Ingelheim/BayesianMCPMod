#' @title simulateData
#' @param n_patients tbd
#'
#' @param dose_levels tbd
#' @param pcov tbd
#' @param SD tbd
#' @param placebo_effect tbd
#' @param max_effect tbd
#' @param mod tbd
#' @param n_simulations tbd
#'
#' @export
simulateData <- function(
    n_patients,
    dose_levels,
    pcov,
    SD,
    placebo_effect,
    max_effect,
    mod,
    n_simulations) {
  ntrt <- length(dose_levels)
  sum_patients <- sum(n_patients)

  simulno <- rep(1:n_simulations,
    each = sum_patients
  )
  ptno <- rep(1:sum_patients,
    times = n_simulations
  )
  dose <- rep(rep(dose_levels, times = n_patients),
    times = n_simulations
  )

  df <- data.frame(simulation = simulno, ptno, dose)
  df <- cbind(
    df,
    DoseFinding::getResp(mod, dose) +
      rep(stats::rnorm(length(dose), mean = 0, sd = SD), times = length(mod))
  )

  return(df)
}
