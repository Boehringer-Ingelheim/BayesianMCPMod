#' @title BayesMCPMod
#'
#' @param ancova tbd
#' @param contMat1 tbd
#' @param crit_prob tbd
#'
#' @return tbd

BayesMCPMod <- function(
    ancova,
    contMat1,
    crit_prob) {
  out <- list()
  contMat <- contMat1[[1]]

  m <- length(contMat[1, ])
  probs <- rep(0, m)
  adjpval <- rep(0, m)
  sign <- 0

  for (j in 1:m) {
    conttheta <- sapply(seq_len(nrow(ancova[["est."]])), function(x) t(contMat[, j]) %*% ancova[["est."]][x, ])
    contvar <- sapply(seq_len(nrow(ancova[["est."]])), function(x) t(contMat[, j]^2) %*% ancova[["cov."]][x, ]^2)
    p_ij <- sapply(seq_len(nrow(ancova[["est."]])), function(x) stats::pnorm(conttheta[x] / sqrt(contvar[x])))


    probs[j] <- sum(ancova[["obs."]] * p_ij)
  }

  max_prob <- max(probs)

  if (max_prob > crit_prob) {
    sign <- 1
  }

  out[["sign."]] <- sign
  out[["adjpval."]] <- max_prob
  out[["posteriorProb"]] <- probs

  return(out)
}

#' @title doAncova
#'
#' @param dose tbd
#' @param response tbd
#' @param prior tbd
#'
#' @return tbd

doAncova <- function(
    dose,
    response,
    prior) {
  mixnorm <- RBesT::mixnorm


  anova <- stats::lm(response ~ factor(dose) - 1)
  mean.sim <- anova$coefficients
  se.mean.sim <- summary(anova)$coefficients[, 2]

  n.dose <- unique(prior[, 1])

  doseAct <- n.dose[-1]

  k <- length(n.dose)

  prior_dose <- lapply(1:k, function(x) matrix(subset(prior, prior[, 1] == n.dose[x])[, -1], ncol = 3))

  comb.prior <- sapply(1:k, function(x) nrow(prior_dose[[x]]))


  post_sep <- list()

  for (i in 1:k) {
    args <- lapply(1:comb.prior[i], function(x) prior_dose[[i]][x, ])
    mixed.prior <- do.call("mixnorm", args)
    post_sep[[i]] <- RBesT::postmix(mixed.prior, m = mean.sim[i], se = se.mean.sim[i])
  }

  comb.post <- prod(comb.prior)

  args <- lapply(1:k, function(x) 1:comb.prior[x])
  comb.ind <- do.call("expand.grid", args)

  post.weight <- matrix(sapply(1:k, function(x) sapply(1:comb.post, function(y) post_sep[[x]][1, comb.ind[y, x]])), nrow = comb.post)
  post.weight <- apply(post.weight, 1, prod)
  post.mean <- matrix(sapply(1:k, function(x) sapply(1:comb.post, function(y) post_sep[[x]][2, comb.ind[y, x]])), nrow = comb.post)
  post.sd <- matrix(sapply(1:k, function(x) sapply(1:comb.post, function(y) post_sep[[x]][3, comb.ind[y, x]])), nrow = comb.post)
  post.sd.overall <- sapply(1:k, function(x) summary(post_sep[[x]])[2])
  post.mean.pbo <- RBesT::qmix(post_sep[[1]], p = 0.5)
  post.mean.act1 <- RBesT::qmix(post_sep[[2]], p = 0.5)
  post.mean.act2 <- RBesT::qmix(post_sep[[3]], p = 0.5)
  post.mean.act3 <- RBesT::qmix(post_sep[[4]], p = 0.5)
  diff <- c(
    post.mean.act1 - post.mean.pbo,
    post.mean.act2 - post.mean.pbo,
    post.mean.act3 - post.mean.pbo
  )
  post.mean.max.diff <- c(max(diff), which.max(diff))

  output <- list()


  output[["obs."]] <- post.weight
  output[["est."]] <- post.mean
  output[["cov."]] <- post.sd
  output[["post.sd."]] <- post.sd.overall
  output[["pbomedian."]] <- post.mean.pbo
  output[["act1median."]] <- post.mean.act1
  output[["act2median."]] <- post.mean.act2
  output[["act3median."]] <- post.mean.act3
  output[["maxdiffmedian+which."]] <- post.mean.max.diff

  return(output)
}