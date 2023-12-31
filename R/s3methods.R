## BayesianMCPMod -----------------------------------------

#' @export
print.BayesianMCPMod <- function (
    
  x,
  ...
  
) {
  
  model_names <- colnames(x$BayesianMCP)[
    grepl("post_probs.", colnames(x$BayesianMCP))] |>
    sub(pattern = "post_probs.", replacement = "", x = _) |>
    gsub(pattern = "\\d", replacement = "", x = _) |>
    unique(x = _)
  
  model_success <- colMeans(do.call(rbind, lapply(x$Mod, function (y) {
    
    if (!is.null(y)) {
      
      sapply(y, function (z) z$significant)
      
    } else {
      
      model_signs        <- rep(FALSE, length(model_names))
      names(model_signs) <- model_names
      
      return (model_signs)
      
    }
    
  })))
  
  print(x$BayesianMCP)
  cat("\n")
  cat("Model Significance Frequencies\n")
  print(model_success, ...)
  
  if (any(!is.na(attr(x$BayesianMCP, "ess_avg")))) {
    
    cat("Average Posterior ESS\n")
    print(attr(x$BayesianMCP, "ess_avg"), ...)
    
  }
  
}

## BayesianMCP --------------------------------------------

#' @export
print.BayesianMCP <- function (
    
  x,
  ...
  
) {
  
  n_sim <- nrow(x)
  
  cat("Bayesian Multiple Comparison Procedure\n")
  
  if (n_sim == 1L) {
    
    attr(x, "crit_prob_adj") <- NULL
    attr(x, "success_rate")  <- NULL
    class(x) <- NULL
    
    print.default(x, ...)
    
  } else {
    
    cat("  Estimated Success Rate: ", attr(x, "successRate"), "\n")
    cat("  N Simulations:          ", n_sim)
    
  }
  
}

## ModelFits ----------------------------------------------

#' @title predict.modelFits
#' @description This function performs model predictions based on the provided
#' model and dose specifications 
#' 
#' @param object A modelFits object containing information about the fitted
#' model coefficients
#' @param doses A vector specifying the doses for which a prediction should be
#' done 
#' @param ... Currently without function
#' @examples
#' posterior_list <- list(Ctrl = RBesT::mixnorm(comp1 = c(w = 1, m = 0, s = 1), sigma = 2),
#'                        DG_1 = RBesT::mixnorm(comp1 = c(w = 1, m = 3, s = 1.2), sigma = 2),
#'                        DG_2 = RBesT::mixnorm(comp1 = c(w = 1, m = 4, s = 1.5), sigma = 2) ,  
#'                        DG_3 = RBesT::mixnorm(comp1 = c(w = 1, m = 6, s = 1.2), sigma = 2) ,
#'                        DG_4 = RBesT::mixnorm(comp1 = c(w = 1, m = 6.5, s = 1.1), sigma = 2))
#' models         <- c("emax", "exponential", "sigEmax", "linear")
#' dose_levels    <- c(0, 1, 2, 4, 8)
#' fit            <- getModelFits(models      = models,
#'                                posterior   = posterior_list,
#'                                dose_levels = dose_levels)
#'                                
#' predict(fit, doses = c(0, 1, 3, 4, 6, 8))
#' 
#' @return a list with the model predictions for the specified models and doses
#' 
#' @export
predict.modelFits <- function (
    
  object,
  doses = NULL,
  ...
  
) {

  predictions <- lapply(object, predictModelFit, doses = doses)
  attr(predictions, "doses") <- doses
  
  return (predictions)
  
}

#' @export
print.modelFits <- function (
    
  x,
  ...
  
) {
  
  n_digits = 1
  
  dose_levels <- x[[1]]$dose_levels
  dose_names  <- names(attr(x, "posterior"))
  
  predictions <- t(sapply(x, function (y) y$pred_values))
  colnames(predictions) <- dose_names
  
  out_table <- data.frame(predictions,
                          mEff = sapply(x, function (y) y$max_effect),
                          gAIC = sapply(x, function (y) y$gAIC),
                          w    = sapply(x, function (y) y$model_weight))
  out_table <- apply(as.matrix(out_table), 2, round, digits = n_digits)
  
  if (!is.null(x[[1]]$significant)) {
    
    out_table <- as.data.frame(out_table)
    out_table$Sign <- sapply(x, function (y) y$significant)
    
  }
  
  model_names <- names(x) |>
    gsub("exponential", "exponential", x = _) |>
    gsub("quadratic",   "quadratic  ", x = _) |>
    gsub("linear",      "linear     ", x = _) |>
    gsub("logistic",    "logistic   ", x = _) |>
    gsub("emax",        "emax       ", x = _) |>
    gsub("sigEmax",     "sigEmax    ", x = _)
  
  # cat("Bayesian MCP Model Fits\n\n")
  cat("Model Coefficients\n")
  for (i in seq_along(model_names)) {
    
    coeff_values <- x[[i]]$coeff
    coeff_names  <- names(coeff_values)
    
    cat(model_names[i],
        paste(coeff_names, round(coeff_values, n_digits), sep = " = "), "\n",
        sep = " ")
    
  }
  cat("\n")
  cat("Dose Levels\n",
      paste(dose_names, round(dose_levels, n_digits), sep = " = "), "\n")
  cat("\n")
  cat("Predictions, Maximum Effect, gAIC, Model Weights & Significance\n")
  print(out_table, ...)
  
}

## plot.ModelFits()
## see file R/plot.R

## postList -----------------------------------------------

#' @export
summary.postList <- function (
    
  object,
  ...
  
) {
  
  summary_list        <- lapply(object, summary, ...)
  names(summary_list) <- names(object)
  summary_tab         <- do.call(rbind, summary_list)
  
  return (summary_tab)
  
}

#' @export
print.postList <- function (
    
  x,
  ...
  
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
  
  summary_tab <- summary.postList(x)
  
  names(x) <- rownames(summary_tab)
  class(x) <- NULL
  
  list_out <- list(summary_tab, getMaxDiff(summary_tab[, 4]), x)
  names(list_out) <- c("Summary of Posterior Distributions",
                       "Maximum Difference to Control and Dose Group",
                       "Posterior Distributions")
  
  print(list_out, ...)
  
}

