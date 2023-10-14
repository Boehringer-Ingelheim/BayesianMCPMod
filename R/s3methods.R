## BayesianMCPMod -----------------------------------------

#' @export
print.BayesianMCPMod <- function (
    
  x,
  ...
  
) {
  
  n_models      <- ncol(x$BayesianMCP) - 2L
  model_names   <- colnames(x$BayesianMCP)[-c(1, 2)] |>
    sub(pattern = "post_probs.", replacement = "", x = _)
  
  model_success <- colMeans(do.call(rbind, lapply(x$Mod, function (y) {
    
    if (!is.null(y)) {
      
      sapply(y, function (z) z$significant)
      
    } else {
      
      model_signs <- rep(FALSE, n_models)
      names(model_signs) <- model_names
      
      return (model_signs)
      
    }
    
  })))
  
  print(x$BayesianMCP)
  cat("\n")
  cat("Model Significance Frequencies\n")
  print(model_success, ...)
  
}

## BayesianMCP --------------------------------------------

#' @export
print.BayesianMCP <- function (
    
  x,
  ...
  
) {
  
  power <- mean(x[, 1])
  n_sim <- nrow(x)
  
  cat("Bayesian Multiple Comparison Procedure\n")
  cat("  Estimated Success Rate: ", power, "\n")
  cat("  N Simulations:          ", n_sim)
  
}

## ModelFits ----------------------------------------------

#' @export
predict.ModelFits <- function (
    
  object,
  doses = NULL,
  ...
  
) {
  
  lapply(object, predictModelFit, doses = doses, ...)
  
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
                          gAIC = sapply(x, function (y) y$gAIC))
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
  cat("Predictions, Maximum Effect, gAIC & Significance\n")
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

