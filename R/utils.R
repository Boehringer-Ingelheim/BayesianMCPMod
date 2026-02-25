

## Make parallel processing optional
useFuture <- function () isTRUE(requireNamespace("future.apply", quietly = TRUE))

optPar_lapply <- function (X, FUN, ...) {
  
  if (useFuture()) {
    
    return (future.apply::future_lapply(
      X, FUN,
      future.seed       = FALSE,
      future.scheduling = 1,
      future.packages   = c("DoseFinding", "RBesT", "nloptr"),
      ...))
    
  } else {
    
    return (lapply(X, FUN, ...))
    
  }
  
}

optPar_apply <- function (X, MARGIN, FUN, ...) {
  
  if (useFuture()) {
    
    return (future.apply::future_apply(
      X, MARGIN, FUN,
      future.seed       = FALSE, 
      future.scheduling = 1,
      future.packages   = c("DoseFinding", "RBesT")),
      ...)
    
  } else {
    
    return (apply(X, MARGIN, FUN, ...))
    
  }
  
}
