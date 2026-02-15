#' @title simulateData
#' 
#' @description
#' Function to simulate patient level data for a normally distributed endpoint
#' @param n_patients Vector containing number of patients as a numerical
#' value per dose-group.
#' @param dose_levels Vector containing the different dosage levels.
#' @param sd Standard deviation on patient level. Can be NULL if `probability_scale` is TRUE. Default NULL.
#' @param mods An object of class "Mods" as specified in the DoseFinding package. Can be NULL if ´dr_means´ is not NULL. Default NULL.
#' @param n_sim Number of simulations to be performed,
#' Default is 1000
#' @param true_model A character for model name, e.g. "emax". Assumed true underlying model.
#' If NULL, all dose-response models included in the mods input parameter will be used. Default NULL.
#' @param dr_means an optional vector, with information about assumed effects per dose group. Default NULL.
#' @param probability_scale A boolean to specify if the trial has a continuous or a binary outcome. Setting to TRUE will transform calculations from the logit scale to the probability scale, which can be desirable for a binary outcome. Default FALSE.
#' 
#' @examples
#' models <- DoseFinding::Mods(linear      = NULL,
#'                             linlog      = NULL,
#'                             emax        = c(0.5, 1.2),
#'                             exponential = 2, 
#'                             doses       = c(0, 0.5, 2,4, 8),
#'                             maxEff      = 6)
#' dose_levels <- c(0, 0.5, 2, 4, 8)
#' sd          <- 12
#' n_patients  <- c(40, 60, 60, 60, 60)
#' 
#' sim_data <- simulateData(n_patients  = n_patients,
#'                          dose_levels = dose_levels,
#'                          sd          = sd,
#'                          mods        = models)
#' 
#' head(sim_data)
#' 
#' # custom response "model" shape
#' custom_dose_response <- c(1, 2, 3, 4, 5)
#' sim_data_custom_dr   <- simulateData(n_patients  = n_patients,
#'                                      dose_levels = dose_levels,
#'                                      sd          = sd,
#'                                      dr_means    = custom_dose_response)
#' 
#' head(sim_data_custom_dr)
#'
#' @return A list object, containing patient level simulated data for all assumed true models.
#' Also providing information about simulation iteration, patient number as well as dosage levels.
#' 
#' @export
simulateData <- function(
    
  n_patients,
  dose_levels,
  sd                = NULL,
  mods              = NULL,
  n_sim             = 1e3,
  true_model        = NULL,
  dr_means          = NULL,
  probability_scale = FALSE
  
) {
  
  checkmate::check_vector(n_patients, any.missing = FALSE, len = length(dose_levels))
  checkmate::check_double(dose_levels, lower = 0, any.missing = FALSE, len = length(n_patients))
  checkmate::check_double(sd, len = 1, null.ok = TRUE, lower = 0, upper = Inf)
  checkmate::check_class(mods, classes = "Mods", null.ok = TRUE)
  checkmate::check_numeric(n_sim, lower = 0, upper = Inf, len = 1)
  checkmate::check_string(true_model, null.ok = TRUE)
  checkmate::assert_flag(probability_scale)
  
  if (!probability_scale) stopifnot("Must provide 'sd' argument for simulation." = !is.null(sd))

  # stopifnot("Either 'dr_means' or 'mods' must be NULL." =
  #             is.null(dr_means) & !is.null(mods) |
  #             !is.null(dr_means) & is.null(mods))
  
  if (!probability_scale) stopifnot(!is.null(sd))
  
  if (!is.null(true_model)) stopifnot(!is.null(mods))
  
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
      data     = rep(unlist(mapply(rep, dr_means, n_patients)), n_sim),
      ncol     = 1L,
      dimnames = list(NULL, "dr_response"))
    
  }
  
  if (probability_scale) {
    
    # simulate responders from binomial distribution
    model_response_rates <- RBesT::inv_logit(model_responses)
    simulated_responses  <- matrix(
      rbinom(length(model_response_rates), size = 1L, prob = as.vector(model_response_rates)),
      nrow = nrow(model_response_rates), ncol = ncol(model_response_rates))
    colnames(simulated_responses) <- colnames(model_response_rates)
    
  } else {
    
    # purposefully simulating random noise per patient, not per patient and model shape
    random_noise        <- stats::rnorm(nrow(sim_info), mean = 0, sd = sd) 
    simulated_responses <- model_responses + random_noise
    
  }
  
  sim_data <- cbind(sim_info, simulated_responses)
  
  rownames(sim_data)                  <- NULL
  attr(sim_data, "probability_scale") <- probability_scale
  
  return (sim_data)
  
}

getModelData <- function (
    
  sim_data,
  model_name
  
) {
  
  model_data <- sim_data[, c("simulation", "dose", model_name)]
  colnames(model_data)[3] <- "response"
  
  attr(model_data, "probability_scale") <- attr(sim_data, "probability_scale")
  
  return (model_data)
  
}
