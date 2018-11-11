#' @title Pseudo inertia of a cluster 
#' @description The pseudo inertia of a cluster is calculated from a dissimilarity matrix
#' and not from a data matrix. 
#' @param D an object of class "dist" with the dissimilarities between the n observations. 
#' The function \code{\link{as.dist}} can be used to transform an object of class matrix to object of class "dist".
#' @param indices a vector with the indices of the subset of observations.
#' @param wt vector with the weights of the n observations
#' @examples 
#' data(estuary)
#' n <- nrow(estuary$dat)
#' Z <- scale(estuary$dat)*sqrt(n/(n-1))
#' inertdiss(dist(Z)) # pseudo inertia
#' inert(Z) #equals for euclidean distance
#' 
#' w <- estuary$map@data$POPULATION # non uniform weights 
#' inertdiss(dist(Z),wt=w)
#' 
#' @references 
#' M. Chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco. ClustGeo: an R package
#' for hierarchical clustering with spatial constraints.
#' Comput Stat (2018) 33: 1799-1822. 
#' 
#' @export
#' 
inertdiss <- function (D, indices=NULL, wt = NULL) {
  n <- as.integer(attr(D, "Size"))
  if (is.null(indices)) indices <- 1:n
  if (is.null(wt)) wt <- rep(1/n,n)
  D <- as.matrix(D)
  if (length(indices) > 1) {
    subD <- D[indices,indices]
    subw <- wt[indices]
    mu <- sum(wt[indices])
    inert <-sweep(subD^2, 1, FUN="*", STATS=subw)
    inert <-sweep(inert, 2, FUN="*", STATS=subw)
    inert <- sum(inert/(2*mu))  
  }
  else {
    inert <- 0
  }
  return(inert)
}

