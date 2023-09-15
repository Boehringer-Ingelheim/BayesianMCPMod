## from DoseFinding
powCalc <- function(alternative, critV, df, corMat, deltaMat, control) {
  mvtnorm.control <- DoseFinding::mvtnorm.control
  
  pmvt <- mvtnorm::pmvt
  
  nC <- nrow(corMat)
  if (alternative[1] == "two.sided") {
    lower <- rep(-critV, nC)
  } else {
    lower <- rep(-Inf, nC)
  }
  upper <- rep(critV, nC)
  if (!missing(control)) {
    if (!is.list(control)) {
      stop("when specified, 'control' must be a list")
    }
    ctrl <- do.call("mvtnorm.control", control)
  } else {
    ctrl <- control
  }
  ctrl$interval <- NULL
  nScen <- ncol(deltaMat)
  res <- numeric(nScen)
  for (i in 1:nScen) {
    pmvtCall <- c(list(lower, upper,
                       df = df, corr = corMat,
                       delta = deltaMat[, i], algorithm = ctrl
    ))
    res[i] <- as.vector(1 - do.call("pmvt", pmvtCall))
  }
  names(res) <- colnames(deltaMat)
  res
}
