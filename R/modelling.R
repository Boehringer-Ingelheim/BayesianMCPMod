#' @title doFit
#' @param dose_levels tbd
#'
#' @param shapes tbd
#' @param posterior tbd
#'
#' @export
doFit <- function(
    dose_levels,
    shapes,
    posterior) {
  mD <- max(dose_levels)
  bounds <- DoseFinding::defBnds(
    mD          = mD,
    emax        = c(0.001, 1.5) * mD,
    exponential = c(0.1, 2) * mD,
    logistic    = matrix(c(0.001, 0.01, 1.5, 1 / 2) * mD, 2)
  )

  summary <- list()

  resp <- c(
    posterior$pbomedian,
    posterior$act1median,
    posterior$act2median,
    posterior$act3median
  )

  cov <- diag(posterior$post.sd)

  for (j in seq_along(shapes)) {
    sublist <- list()
    shape <- shapes[j]
    fit <- DoseFinding::fitMod(
      dose = dose_levels,
      resp = resp,
      S = cov,
      model = shape,
      type = "general",
      bnds = bounds[[shape]]
    )
    fit_vec <- stats::predict(fit, predType = "ls-means")
    gAIC <- DoseFinding::gAIC(fit)
    max_effect <- max(fit_vec) - min(fit_vec)

    sublist[["model"]] <- shape
    sublist[["fit"]] <- fit$coefs
    sublist[["pred"]] <- fit_vec
    sublist[["gAIC"]] <- gAIC
    sublist[["max_effect"]] <- max_effect

    summary[[j]] <- sublist
  }

  return(summary)
}

#' @title estimateModels
#' @param models tbd
#'
#' @param post tbd
#' @param doses tbd
#'
#' @export
estimateModels <- function(
    models,
    post,
    doses) {
  fitModel <- function(model, post, doses) {
    getF <- function(
        theta,
        post,
        dose,
        expression_i) {
      L <- length(post$obs)
      comp <- rep(0, L)

      for (i in 1:L) {
        comp[i] <- eval(expression_i)
      }

      return(sum(post$obs * comp))
    }

    switch(model,
      "emax" = {
        x0 <- c(0, 1, 0.5)
        maxeval <- 5000
        lb <- c(-Inf, -Inf, 0.001 * max(doses))
        ub <- c(Inf, Inf, 1.5 * max(doses))
        expression_i <- quote(sum((post$est[i, ] - (theta[1] + (theta[2] * dose) / (theta[3] + dose)))^2 / (post$cov[i, ]^2)))
      },
      "sigEmax" = {
        x0 <- c(0, 1, 0.5, 1)
        maxeval <- 10000
        lb <- c(-Inf, -Inf, 0.001 * max(doses), 0.5)
        ub <- c(Inf, Inf, 1.5 * max(doses), 10)
        expression_i <- quote(sum((post$Mean[i, ] - (theta[1] + (theta[2] * dose^theta[4]) / (theta[3]^theta[4] + dose^theta[4])))^2 / (post$cov[i, ]^2)))
      },
      "exponential" = {
        x0 <- c(-1, 1, 20)
        maxeval <- 10000
        lb <- c(-Inf, -Inf, 0.1 * max(doses))
        ub <- c(Inf, Inf, 2 * max(doses))
        expression_i <- quote(sum((post$est[i, ] - (theta[1] + theta[2] * (exp(dose / theta[3]) - 1)))^2 / (post$cov[i, ]^2)))
      },
      "quadratic" = {
        x0 <- c(0, 1, -1)
        maxeval <- 10000
        lb <- NULL
        ub <- c(Inf, Inf, 0)
        expression_i <- quote(sum((post$est[i, ] - (theta[1] + theta[2] * dose + theta[3] * dose^2))^2 / (post$cov[i, ]^2)))
      },
      "linear" = {
        x0 <- c(0, 1)
        maxeval <- 10000
        lb <- NULL
        ub <- NULL
        expression_i <- quote(sum((post$est[i, ] - (theta[1] + theta[2] * dose))^2 / (post$cov[i, ]^2)))
      },
      "logistic" = {
        x0 <- c(0, 1, 0.5, 5)
        maxeval <- 10000
        lb <- c(-Inf, -Inf, 0.001 * max(doses), 0.01 * max(doses))
        ub <- c(Inf, Inf, 1.5 * max(doses), 0.5 * max(doses))
        expression_i <- quote(sum((post$est[i, ] - (theta[1] + theta[2] / (1 + exp((theta[3] - dose) / theta[4]))))^2 / (post$cov[i, ]^2)))
      },
      {
        stop("model must be one of 'emax', 'sigEmax', 'exponential', 'linear', 'logistic', 'quadratic'")
      }
    )
    res <- nloptr::nloptr(
      x0 = x0,
      eval_f = function(theta) {
        getF(theta,
          post         = post,
          dose         = doses,
          expression_i = expression_i
        )
      },
      opts = list(
        "algorithm" = "NLOPT_LN_NELDERMEAD",
        maxeval = maxeval
      ),
      lb = lb,
      ub = ub
    )

    return(list(
      model = model,
      fit = res
    ))
  }

  results <- lapply(models, function(model) {
    fitModel(model, post, doses)
  })

  return(results)
}

#' @title getGenAICs
#' @param post tbd
#'
#' @param fit_mods tbd
#' @param doses tbd
#'
#' @export
getGenAICs <- function(
    post,
    fit_mods,
    doses) {
  getGenAIC <- function(
      post,
      fit_mod,
      doses) {
    model <- fit_mod$model
    p <- length(fit_mod$fit$solution)
    theta <- fit_mod$fit$solution
    L <- length(post$obs)

    fitted <- switch(model,
      "emax" = {
        DoseFinding::emax(doses, theta[1], theta[2], theta[3])
      },
      "sigEmax" = {
        DoseFinding::sigEmax(doses, theta[1], theta[2], theta[3], theta[4])
      },
      "exponential" = {
        DoseFinding::exponential(doses, theta[1], theta[2], theta[3])
      },
      "quadratic" = {
        DoseFinding::quadratic(doses, theta[1], theta[2], theta[3])
      },
      "linear" = {
        DoseFinding::linear(doses, theta[1], theta[2])
      },
      "logistic" = {
        DoseFinding::logistic(doses, theta[1], theta[2], theta[3], theta[4])
      },
      {
        stop("model must be one of 'emax', 'sigEmax', 'exponential', 'linear', 'logistic', 'quadratic'")
      }
    )

    gAIC <- rep(0, L)

    for (i in 1:L) {
      gAIC[i] <- sum(1 / post$cov[i, ]^2 * (post$est[i, ] - fitted)^2) + 2 * p
    }

    avgAIC <- stats::weighted.mean(gAIC, w = post$obs)

    res <- list(
      model = model,
      fitted = fitted,
      gAIC = gAIC,
      avgAIC = avgAIC,
      weightedAIC = NULL
    )

    return(res)
  }

  results <- lapply(fit_mods, function(fit_mod) {
    getGenAIC(post, fit_mod, doses)
  })

  avgAIC_values <- sapply(results, function(res) res$avgAIC)
  exp_values <- exp(-0.5 * avgAIC_values)
  denominator <- sum(exp_values)

  for (i in 1:length(results)) {
    results[[i]]$weightedAIC <- exp_values[i] / denominator
  }

  return(results)
}
