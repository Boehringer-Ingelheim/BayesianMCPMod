#' @title getPosterior
#' 
#' @param data tbd
#' @param prior prior_list
#'
#' @export
getPosterior <- function(
    
  data,
  prior_list,
  mu_hat = NULL,
  sd_hat = NULL
  
) {
  
  lapply(split(data, data$simulation), getPosteriorI,
         prior_list = prior_list,
         mu_hat     = mu_hat,
         sd_hat     = sd_hat)
  
}

getPosteriorI <- function(
    
  data_i,
  prior_list,
  mu_hat = NULL,
  sd_hat = NULL
  
) {
  
  if (is.null(mu_hat) && is.null(sd_hat)) {
    
    anova_res <- lm(data_i$response ~ factor(data_i$dose) - 1)
    mu_hat    <- summary(anova_res)$coefficients[, 1]
    sd_hat    <- summary(anova_res)$coefficients[, 2]
    
  } else if (!is.null(mu_hat) && !is.null(sd_hat)) {
    
    stopifnot("m_hat length must match number of dose levels" = 
                length(prior_list) == length(mu_hat),
              "sd_hat length must match number of dose levels" = 
                length(prior_list) == nrow(sd_hat))
    
  } else {
    
    stop ("Both mu_hat and S_hat must be provided.")
    
  }
  
  post_list <- mapply(RBesT::postmix, prior_list, m = mu_hat, se = sd_hat)
  
  names(post_list) <- c("Ctr", paste0("DG_", seq_along(post_list[-1])))
  class(post_list) <- "postList"
  
  return (post_list)
  
}

getPostCombsI <- function (
    
  posterior_i
  
) {
  
  post_params <- list(
    weights = lapply(posterior_i, function (x) x[1, ]),
    means   = lapply(posterior_i, function (x) x[2, ]),
    vars    = lapply(posterior_i, function (x) x[3, ]^2))
  
  post_combs         <- lapply(post_params, expand.grid)
  post_combs$weights <- apply(post_combs$weights, 1, prod)
  
  return (post_combs)
  
}

summary.postList <- function (
    
  post_list,
  ...
  
) {
  
  summary_list        <- lapply(post_list, summary, ...)
  names(summary_list) <- names(post_list)
  summary_tab         <- do.call(rbind, summary_list)
  
  return (summary_tab)
  
}

print.postList <- function (
    
  post_list
  
) {
  
  getMaxDiff <- function (
    
    medians
    
  ) {
    
    diffs <- medians - medians[1]
    
    max_diff       <- max(diffs)
    max_diff_level <- which.max(diffs) - 1
    
    out <- c(max_diff, max_diff_level)
    names(out) <- c("max_diff", "DG")
    
    return (out)
    
  }
  
  summary_tab <- summary.postList(post_list)
  
  names(post_list) <- rownames(summary_tab)
  class(post_list) <- NULL
  
  list_out <- list(summary_tab, getMaxDiff(summary_tab[, 4]), post_list)
  names(list_out) <- c("Summary of Posterior Distributions",
                       "Maximum Difference to Control and Dose Group",
                       "Posterior Distributions")
  
  print(list_out)
  
}
