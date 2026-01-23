padStrings <- function(strings) {
  
  max_length     <- max(nchar(strings))
  padded_strings <- sapply(strings, function(x) {
    
    paste0(x, strrep(" ", max_length - nchar(x)))
    
  })
  
  return (padded_strings)
  
}

printMatrixWithPrefix <- function (
    
  mat,
  prefix = "  "
  
) {
  
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  
  if (!is.null(row_names)) {
    
    if (!is.null(col_names)) {
      
      cat(rep(" ", max(nchar(row_names))), " ", sep = "")
      
    }
    
    row_names <- padStrings(row_names)
    
  }
  
  print_width <- max(nchar(round(mat, options()$digits)),
                     nchar(col_names), na.rm = TRUE)
  
  if (!is.null(col_names)) {
    
    cat(prefix)
    cat(format(col_names, width = print_width, justify = "right"))
    cat("\n")
    
  }
  
  for (i in seq_len(nrow(mat))) {
    
    cat(prefix)
    
    if (!is.null(row_names)) {
      
      cat(row_names[i], " ", sep = "")
      
    }
    
    cat(format(mat[i, ], width = print_width), "\n")
    
  }
  
}

shortenModelNames <- function (model_names, pad_string = FALSE) {
  
  model_names_out <- model_names |>
    gsub("exponential", "exp", x = _) |>
    gsub("quadratic",  "quad", x = _) |>
    gsub("linear",      "lin", x = _) |>
    gsub("logistic",    "log", x = _) |>
    gsub("sigEmax",    "sigE", x = _) |>
    gsub("betaMod",   "betaM", x = _) |>
    gsub("quadratic",  "quad", x = _)
  
  if (pad_string) {
    
    model_names_out <- padStrings(model_names_out)
    
  }
  
  return (model_names_out)
  
}

## BayesianMCPMod -----------------------------------------

#' @export
print.BayesianMCPMod <- function (
    
  x,
  n_mods = 0L,
  ...
  
) {
  
  print(x$BayesianMCP, ...)
  
  if (any(!is.na(attr(x, "MED")))) {
    
    sign_indx <- as.logical(x$BayesianMCP[, "sign"])
    med_info  <- attr(x, "MED")
    
    cat("MED Assessment")
    
    if (any(sign_indx)) {
      
      probability_scale <- attr(x$Mod[[which(sign_indx)[1]]], "probability_scale")
      
      if (probability_scale) cat(" on Probability Scale")
      
    }
    
    cat("\n")
    
    cat("  Selection Method:   ", attr(x, "MEDSelection"), "\n")
    cat("  Identification Rate:", mean(med_info[, "med_reached"]), "\n")
    
    if (!any(sign_indx)) {
      
      cat("    No model shapes are significant\n")
      
    } else {
      
      # Frequency of Test not passed
      freq_not_passed <- 1 - mean(sign_indx)
      
      # Frequency of MED not reached in passed tests
      freq_not_reached <- mean(as.numeric(!med_info[, "med_reached"]) * as.numeric(sign_indx))
      
      # MED Dose Frequencies
      dose_levels   <- x$Mod[[which(sign_indx)[1]]][[1]]$dose_levels[-1]
      med_table     <- table(factor(med_info[, "med"], levels = dose_levels), useNA = "no")
      
      med_mat <- matrix(data  = c(as.numeric(names(med_table)),
                                  as.vector(med_table) / nrow(med_info)),
                        nrow  = 2,
                        byrow = TRUE)
      rownames(med_mat) <- c("Dose Level:", "MED Freq:")
      
      printMatrixWithPrefix(med_mat, prefix = "   ")
      
      cat("  MED not reached Freq:       ", freq_not_reached, "\n")
      cat("  No success in MCP step Freq:", freq_not_passed, "\n")
      
    }
    
  } 
  
  if (length(x$Mod) == 1) {

    print (x$Mod[[1]], ...)

  } else {
    
    for (m in seq_len(n_mods)) {
      
      print (x$Mod[[m]], ...)
      
    }
    
  }
  
  invisible(x)
  
}

## BayesianMCP --------------------------------------------

#' @export
print.BayesianMCP <- function(x, ...) {
  
  n_sim <- nrow(x)
  
  cat("Bayesian Multiple Comparison Procedure\n")
  
  if (n_sim == 1L) {
    
    cat("  Significant:                  ", x[1, "sign"], "\n")
    cat("  Critical Probability:         ", x[1, "crit_prob_adj"], "\n")
    cat("  Maximum Posterior Probability:", x[1, "max_post_prob"], "\n")
    
    cat("Posterior Probabilities for Model Shapes\n")
    
    model_probs <- x[1, grep("^post_probs\\.", colnames(x))]
    model_mat   <- matrix(
      data     = model_probs,
      nrow     = 1,
      dimnames = list(
        "Posterior Prob",
        shortenModelNames(gsub(
          pattern     = "post_probs\\.",
          replacement = "",
          x           = names(model_probs)))))
    
    printMatrixWithPrefix(
      rbind(model_mat, Significant = x[, -c(1, 2, 3)] > attr(x, "critProbAdj")))
    
    if (any(!is.na(attr(x, "essAvg")))) {
      
      cat("Average Posterior ESS\n")
      
      ess_vec       <- as.vector(attr(x, "essAvg"))
      ess_vec_names <- names(attr(x, "essAvg"))
      
      print_width <- max(nchar(ess_vec), nchar(ess_vec_names), na.rm = TRUE)
      
      cat("  Dose Level:  ", format(ess_vec_names, width = print_width, justify = "right"), "\n")
      cat("  Avg Post ESS:", format(ess_vec, width = print_width), "\n")
      
    }
    
  } else {
    
    model_successes <- getModelSuccesses(x)
    
    cat("  Estimated Success Rate:", attr(x, "successRate"), "\n")
    cat("  N Simulations:         ", n_sim, "\n")
    
    m_succ       <- as.vector(model_successes)
    m_succ_names <- shortenModelNames(names(model_successes))
    
    print_width  <- max(nchar(m_succ), nchar(m_succ_names), na.rm = TRUE)
    
    cat("   Model Shape:      ", format(m_succ_names, width = print_width, justify = "right"), "\n")
    cat("   Significance Freq:", format(m_succ, width = print_width), "\n")
    
  }
  
  invisible(x)
  
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
#' models         <- c("emax", "exponential", "sigEmax", "linear", "betaMod")
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
  doses             = NULL,
  probability_scale = attr(object, "probability_scale"),
  ...
  
) {
  
  checkmate::assert_double(doses, lower = 0, any.missing = FALSE, unique = TRUE,
                           sorted = TRUE, finite = TRUE, null.ok = TRUE)
  checkmate::assert_flag(probability_scale)
  
  model_fits  <- object
  model_names <- names(model_fits)
  
  warning_displayed <- FALSE
  withCallingHandlers(
  
    expr = {
      
      predictions <- lapply(model_fits[model_names != "avgFit"],
                            predictModelFit, doses = doses)
      
      if ("avgFit" %in% model_names) {
        
        preds_avg_fit <- predictAvgFit(model_fits, doses = doses) # predictAvgFit calls predict.modelFits !!
        predictions   <- c(list(avgFit = preds_avg_fit), predictions)
        
      }
      
    },
    
    warning = function(w) {
      
      if (grepl("The specified range exceeds the bounds of the original dose range.", conditionMessage(w))) {
        
        if (warning_displayed) {
          
          invokeRestart("muffleWarning")
          
        } else {
          
          warning_displayed <<- TRUE
          
        }
        
      }
      
    }
      
  )
  
  if (probability_scale) {
    
    predictions <- lapply(predictions, RBesT::inv_logit)
    
  }
  
  attr(predictions, "doses") <- doses
  
  return (predictions)
  
}

#' @export
print.modelFits <- function (
    
  x,
  n_digits          = 1,
  probability_scale = attr(x, "probability_scale"),
  ...
  
) {
  
  checkmate::assert_integerish(n_digits, lower = 0)
  checkmate::assert_flag(probability_scale)
  
  if (probability_scale & n_digits < 2) n_digits <- 2
  
  ## Model Coefficients
  
  model_names <- shortenModelNames(names(x), pad_string = TRUE)
  
  cat("Model Coefficients\n")
  
  for (i in seq_along(x)) {
    
    if (x[[i]]$model != "avgFit") {
      
      coeff_values <- x[[i]]$coeff
      coeff_names  <- names(coeff_values)
      
      cat(" ", model_names[i],
          paste(coeff_names, round(coeff_values, n_digits),
                sep      = " = ",
                collapse = ", "), "\n",
          sep = " ")
      
    }
    
  }
  
  ## Dose Levels
  
  dose_levels <- x[[1]]$dose_levels
  dose_names  <- names(attr(x, "posterior"))
  
  cat("Dose Levels\n", "", # "" is required for proper indentation
      paste(dose_names, round(dose_levels, n_digits),
            sep      = " = ",
            collapse = ", "), "\n")
  
  ## Predictions, Maximum Effect, gAIC & avgFit Model Weights
  
  cat("Predictions")
  
  if (probability_scale) {
    
    cat(" on Prob. Scale")
    
  }
  
  cat(", Maximum Effect, gAIC")
  
  predictions           <- t(sapply(x, function (y) y$pred_values))
  colnames(predictions) <- dose_names
  
  out_table <- data.frame(predictions,
                          mEff = sapply(x, function (y) y$max_effect),
                          gAIC = sapply(x, function (y) y$gAIC),
                          w    = sapply(x, function (y) y$model_weight))
  
  if (!is.null(x[[1]]$significant)) {
    
    model_sig      <- TRUE
    out_table$Sign <- sapply(x, function (y) y$significant)
    
  } else {
    
    model_sig <- FALSE
    
  }
  
  if (model_sig) {
    
    cat(", avgFit Model Weights & Significance\n")
    
  } else {
    
    cat(" & avgFit Model Weights\n")
    
  }
  
  out_table <- apply(as.matrix(out_table), 2, round, digits = n_digits)
  rownames(out_table) <- shortenModelNames(rownames(out_table))
  
  printMatrixWithPrefix(out_table)
  
  invisible(x)
  
}

## plot.ModelFits()
## see file R/plot.R

## postList -----------------------------------------------

#' @export
summary.postList <- function (
    
  object,
  probability_scale = FALSE,
  ...
  
) {
  
  checkmate::assert_flag(probability_scale)
  
  summary_list        <- lapply(object, summary, ...)
  names(summary_list) <- names(object)
  summary_tab         <- do.call(rbind, summary_list)
  
  ## Transformations for probability scale, X ~ N(mu, sigma^2), p ~ expit(X) = inv_logit(X)
  if (probability_scale) {
    
    ## quantiles - simply use expit(), i.e., inv_logit()
    q_indices                <- setdiff(colnames(summary_tab), c("mean", "sd"))
    summary_tab[, q_indices] <- apply(summary_tab[, q_indices], 2, RBesT::inv_logit)
    
    ## means and standard deviations - numerical integration
    getNormal2LogisticNormalMoments <- function (mu, sigma) {
      
      # integrand for the mean E[p]
      integrand_mean <- function (x) {
        
        p <- RBesT::inv_logit(x)
        
        return (p * dnorm(x, mean = mu, sd = sigma))
        
      }
      
      # integrand for second moment E[p^2]
      integrand_second <- function (x) {
        
        p <- RBesT::inv_logit(x)
        
        return (p^2 * dnorm(x, mean = mu, sd = sigma))
        
      }
      
      # numerical integration
      m1 <- integrate(integrand_mean,   -Inf, Inf, rel.tol = 1e-9)$value
      m2 <- integrate(integrand_second, -Inf, Inf, rel.tol = 1e-9)$value
      
      # standard deviation
      sd_p <- sqrt(m2 - m1^2)
      
      return (c(mean = m1, sd = sd_p))
      
    }
    
    summary_tab[, c("mean", "sd")] <- t(mapply(getNormal2LogisticNormalMoments,
                                               mu    = summary_tab[, "mean"],
                                               sigma = summary_tab[, "sd"]))
    
  }
  
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
  
  ## removes the covariance matrices from the print out
  attr(x, "posteriorInfo") <- NULL
  
  list_out <- list(summary_tab, getMaxDiff(summary_tab[, 4]), x)
  
  names(list_out) <- c("Summary of Posterior Distributions",
                       "Maximum Difference to Control and Dose Group",
                       "Posterior Distributions")
  
  print(list_out, ...)
  cat("Note: For a more detailed posterior output including\n",
      "     covariance matricies across mixture components,\n",
      "     please use attr(x, 'posteriorInfo').")
  
}
