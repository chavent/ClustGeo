#' @title Dissimilarity based pseudo within-cluster inertia of a partition
#' @description This function performs the pseudo within-cluster inertia of a partition from a dissimilarity matrix. 
#' @param D an object of class "dist" with the dissimilarities between the n observations. 
#' The function \code{\link{as.dist}} can be used to transform an object of class matrix to object of class "dist".
#' @param part a vector with group membership.
#' @param wt vector with the weights of the observations
#' @references 
#' M.chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco.  ClustGeo:  an R package 
#' for hierarchical clustering with spatial constraints	arXiv:1707.03897 [stat.CO]
#' @export
withindiss <- function (D, part, wt = NULL) {
  n <- as.integer(attr(D, "Size"))
  if (is.null(wt)) wt <- rep(1/n,n)
  k <-length(unique(part))
  W <- 0
  for (i in 1:k)
  {
    A <- which(part==i)
    W <- W + inertdiss(D,A,wt)
  }
  return(W)
}