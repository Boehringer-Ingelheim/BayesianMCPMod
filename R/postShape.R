#' @title postShape
#' @param data tbd
#'
#' @param prior tbd
#' @param n_simulations tbd
#'
#' @export
postShape <- function(
    data,
    prior,
    n_simulations) {
  ancova <- list()
  for (i in seq_len(n_simulations)) {
    datai <- data[data$simulation == i, ]
    dosei <- datai$dose
    responsei <- datai$response
    ancovai <- doAncova(
      dose = dosei,
      response = responsei,
      prior = prior
    )
    ancova[[i]] <- ancovai
  }

  return(ancova)
}
