getMCMCSample <- function (

  model_fits,
  n_samples = 1e4

) {

  post_list      <- attr(model_fits, "posterior")
  
  mu_hat_samples <- sapply(post_list, RBesT::rmix, n = n_samples)
  sd_hat         <- summary.postList(post_observed)[, 2]^2
  
  getModelFitOpt()

  return ()

}