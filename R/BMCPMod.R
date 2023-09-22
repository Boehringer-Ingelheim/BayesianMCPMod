#' @title BMCPMod
#' 
#' @param posteriors_list tbd
#' @param contr_mat tbd
#' @param crit_prob tbd
#' 
#' @export
BMCPMod <- function(
    
  posteriors_list,
  contr_mat,
  crit_prob
  
) {
  
  lapply(posteriors_list, BayesMCPMod, contr_mat, crit_prob)
  
}

BayesMCPMod <- function (
    
  posterior_i,
  contr_mat,
  crit_prob
  
) {
  
  getPostProb <- function (
    
    contr_j,     # j: dose level
    post_combs_i # i: simulation outcome
    
  ) {
    
    ## Test statistic = sum over all components of
    ## posterior weight * normal probability distribution of
    ## critical values for doses * estimated mean / sqrt(product of critical values for doses)
    
    ## Calculation for each component of the posterior
    contr_theta   <- apply(post_combs_i$means, 1, `%*%`, contr_j)
    contr_var     <- apply(post_combs_i$vars, 1, `%*%`, contr_j^2)
    contr_weights <- post_combs_i$weights
    
    ## P(c_m * theta > 0 | Y = y) for a shape m (and dose j)
    post_probs <- sum(contr_weights * stats::pnorm(contr_theta / sqrt(contr_var)))
    
    return (post_probs)
    
  }
  
  post_combs_i <- getPostCombsI(posterior_i)
  post_probs   <- apply(contr_mat$contMat, 2, getPostProb, post_combs_i)
  
  res_list <- list(
    sign       = ifelse(max(post_probs) > crit_prob, 1, 0),
    p_val      = max(post_probs),
    post_probs = post_probs)
  
  return (res_list)
  
}
