#' @title Inertia of a cluster
#' @description Computes the inertia of a cluster i.e. on a subset of rows of a data matrix.
#' @param Z matrix data
#' @param indices vectors representing the subset of rows
#' @param wt weight vector
#' @param M diagonal distance matrix
#' 
#' @examples 
#' data(estuary)
#' n <- nrow(estuary$dat)
#' Z <- scale(estuary$dat)*sqrt(n/(n-1))
#' inert(Z) # number of variables
#' 
#' w <- estuary$map@data$POPULATION # non uniform weights 
#' inert(Z,wt=w)
#' 
#' @export

inert <- function (Z, indices=1:nrow(Z), wt = rep(1/nrow(Z), nrow(Z)), M = rep(1, ncol(Z))) {
  # checks there is not a empty class
  if (is.null(indices)) stop('inert : the cluster is empty!')
  if (length(indices) > 1) {
    subZ <- Z[indices,]
    mu <- sum(wt[indices])
    g <- colSums(wt[indices] * Z[indices,]) / mu
    sqdistg <- function(x) {sum(M*(x-g)^2)}
    inert <- sum(wt[indices]*apply(subZ,1,sqdistg))
  }
  else {
    inert <- 0
  }
  return(inert)
}