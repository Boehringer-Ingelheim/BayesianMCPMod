plot.estMod <- function (
    
  est_mod
  # dose_levels
  # posteriors = posterior_linear[[2]]
  
) {
  
  model <- est_mod$model
  theta <- est_mod$fit$solution
  
  dose  <- seq(min(dose_levels), max(dose_levels), length.out = 1e3)
  
  switch(model,
         "emax" = {
           resp_expr <- quote(theta[1] + (theta[2] * dose) / (theta[3] + dose))},
         "sigEmax" = {
           resp_expr <- quote(theta[1] + (theta[2] * dose^theta[4]) / (theta[3]^theta[4] + dose^theta[4]))},
         "exponential" = {
           resp_expr <- quote(theta[1] + theta[2] * (exp(dose / theta[3]) - 1))},
         "quadratic" = {
           resp_expr <- quote(theta[1] + theta[2] * dose + theta[3] * dose^2)},
         "linear" = {
           resp_expr  <- quote(theta[1] + theta[2] * dose)},
         "logistic" = {
           resp_expr <- quote(theta[1] + theta[2] / (1 + exp((theta[3] - dose) / theta[4])))},
         {
           stop(GENERAL$ERROR$MODEL_OPTIONS)}
  )
  
  df <- data.frame(dose = dose, response = eval(resp_expr))
  
  data.frame(dose = dose_levels, obs = )
  
  plt <- ggplot2::ggplot(data = df) +
    ggplot2::geom_line(ggplot2::aes(dose, response)) +
    ggplot2::geom_point()
  
  return(plt)
  
}

plot.estMods <- function (
    
  est_mods
  
) {
  
  plts <- lapply(est_mods, plot)
  
  
  
}
