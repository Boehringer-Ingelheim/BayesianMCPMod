#' @title simulateData
#' 
#' @param n_patients tbd
#' @param dose_levels tbd
#' @param sd tbd
#' @param mods tbd
#' @param n_sim tbd
#' @param true_model tbd
#' 
#' @export
simulateData <- function(
    
  n_patients,
  dose_levels,
  sd,
  mods,
  n_sim = 1e3,
  true_model = NULL
  
) {
  
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
  
  model_responses <- DoseFinding::getResp(mods, sim_info$dose)
  random_noise    <- stats::rnorm(nrow(sim_info), mean = 0, sd = sd)
  
  sim_data <- cbind(sim_info, model_responses + random_noise)
  
  if (!is.null(true_model)) {
    
    sim_data <- getModelData(sim_data, true_model)
    
  }
  
  return (sim_data)
  
}

getModelData <- function (
    
  sim_data,
  model_name
  
) {
  
  model_data              <- sim_data[, c("simulation", "dose", model_name)]
  colnames(model_data)[3] <- "response"
  
  return (model_data)
  
}
