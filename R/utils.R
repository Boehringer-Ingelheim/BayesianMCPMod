

## Make parallel processing optional
useFuture <- function () isTRUE(requireNamespace("future.apply", quietly = TRUE))

optPar_lapply <- function (X, FUN, ...) {
  
  if (useFuture()) {
    
    return (future.apply::future_lapply(X, FUN, ...))
    
  } else {
    
    return (lapply(X, FUN, ...))
    
  }
  
}

optPar_apply <- function (X, MARGIN, FUN, ...) {
  
  if (useFuture()) {
    
    return (future.apply::future_apply(X, MARGIN, FUN, ...))
    
  } else {
    
    return (apply(X, MARGIN, FUN, ...))
    
  }
  
}
