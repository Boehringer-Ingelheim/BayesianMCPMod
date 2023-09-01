#' @title BMCPMod
#' @param ancova1 tbd
#'
#' @param cont_Mat1 tbd
#' @param crit_prob tbd
#' @param n_simulations tbd
#'
#' @export
BMCPMod <- function(
    ancova1,
    cont_Mat1,
    crit_prob,
    n_simulations) {
  adj1_p <- list()

  for (i in 1:n_simulations) {
    ancova <- ancova1[[i]]
    adj1_p[[i]] <- BayesMCPMod(
      ancova,
      cont_Mat1,
      crit_prob
    )
  }

  return(adj1_p)
}
