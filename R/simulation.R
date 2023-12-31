#' @title simulateData
#' 
#' @description
#' Function to simulate patient level data for a normally distributed endpoint
#' @param n_patients Vector containing number of patients as a numerical
#' value per dose-group.
#' @param dose_levels Vector containing the different dosage levels.
#' @param sd Standard deviation on patient level.
#' @param mods An object of class "Mods" as specified in the DoseFinding package.
#' @param n_sim Number of simulations to be performed,
#' Default is 1000
#' @param true_model Default value is NULL.
#' Assumed true underlying model. Provided via a String. e.g. "emax".
#' In case of NULL, all dose-response models, included in the mods input parameter will be used.
#' @param dr_means a vector, with information about assumed effects per dose group. Default NULL.
#' 
#' @examples
#' models <- DoseFinding::Mods(linear      = NULL,
#'                             linlog      = NULL,
#'                             emax        = c(0.5, 1.2),
#'                             exponential = 2, 
#'                             doses       = c(0, 0.5, 2,4, 8),
#'                             maxEff      = 6)
#' dose_levels <- c(0, 0.5, 2,4, 8)
#' sd          <- 12
#' n_patients  <- c(40, 60, 60, 60, 60)
#' 
#' sim_data <- simulateData(n_patients  = n_patients,
#'                          dose_levels = dose_levels,
#'                          sd          = sd,
#'                          mods        = models,
#'                          n_sim       = 100)
#' 
#' sim_data
#'
#' @return A list object, containing patient level simulated data for all assumed true models.
#' Also providing information about simulation iteration, patient number as well as dosage levels.
#' 
#' @export
simulateData <- function(
    
  n_patients,
  dose_levels,
  sd,
  mods,
  n_sim      = 1e3,
  true_model = NULL,
  dr_means   = NULL
  
) {
  
  checkmate::check_vector(n_patients, any.missing = FALSE, len = length(dose_levels))
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(n_patients))
  checkmate::check_double(sd, len = 1, null.ok = FALSE, lower = 0, upper = Inf)
  checkmate::check_class(mods, classes = "Mods")
  checkmate::check_numeric(n_sim, lower = 0, upper = Inf, len = 1)
  checkmate::check_string(true_model, null.ok = TRUE)
  
  if (!is.null(true_model)) {
    
    n_sim            <- 1
    mods_attr        <- attributes(mods)
    mods_attr$names  <- true_model
    mods             <- mods[true_model]
    attributes(mods) <- mods_attr
    
  }
  
  sim_info <- data.frame(
    simulation = rep(seq_len(n_sim), each = sum(n_patients)),
    ptno       = rep(seq_len(sum(n_patients)), times = n_sim),
    dose       = rep(rep(dose_levels, times = n_patients), times = n_sim))
  
  if (is.null(dr_means)) {
    
    model_responses <- DoseFinding::getResp(mods, sim_info$dose)
    
  } else {
    
    stopifnot(identical(length(dose_levels), length(dr_means)))
    
    model_responses <- matrix(
      data     = rep(unlist(mapply(rep, dr_means, n_patients)),
                     n_sim * length(mods)),
      ncol     = length(mods),
      dimnames = list(NULL, names(mods)))
    
  }
  
  random_noise <- stats::rnorm(nrow(sim_info), mean = 0, sd = sd)
  sim_data     <- cbind(sim_info, model_responses + random_noise)
  
  if (!is.null(true_model)) {
    
    sim_data <- getModelData(sim_data, true_model)
    
  }
  
  if (!is.null(dr_means)) {
    
    sim_data           <- sim_data[1:4]
    colnames(sim_data) <- c(colnames(sim_data)[-4], "dose_resp")
    
  }
  
  return (sim_data)
  
}

getModelData <- function (
    
  sim_data,
  model_name
  
) {
  
  model_data <- sim_data[, c("simulation", "dose", model_name)]
  colnames(model_data)[3] <- "response"
  
  return (model_data)
  
}
