#' @title Ward aggregation measures between singletons
#' @description This function calculates the Ward aggregation measures between pairs of singletons. 
#' @param D a object of class "dist" with the dissimilarities between the n obsevations. 
#' The function \code{\link{as.dist}} can be used to transform an object of class matrix to object of class "dist".
#' @param wt vector with the weights of the observations. By default, wt=NULL corresponds to
#' the case where all observations are weighted by 1/n.
#' @return Returns an object of class dist with the Ward aggregation measures between the n singletons.
#' @details The Ward agreggation measure between to singletons i and j weighted by wi and wj is : (wiwj)/(wi+wj)dij^2
#' where dij is the dissimilarity between i and j. 
#' @references 
#' M.chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco.  ClustGeo:  an R package 
#' for hierarchical clustering with spatial constraints	arXiv:1707.03897 [stat.CO]
#' @export

wardinit <- function(D,wt=NULL) {
  n <- as.integer(attr(D, "Size"))
  if (is.null(n)) 
    stop("invalid dissimilarities",call.=FALSE)
  if (is.null(wt)) delta <- D^2/(2*n) else
  {
    delta <-  as.matrix(D)
    delta <-sweep(delta^2, 1, FUN="*", STATS=wt)
    delta <-sweep(delta, 2, FUN="*", STATS=wt)
    S <- matrix(rep(wt,n),n,n)
    delta <- delta/(S+t(S))
  }
  return(as.dist(delta))
} 

# ward_init <- function(d,wt) {
#   n <- as.integer(attr(d, "Size"))
#   diss <-  d               
#   if (is.null(wt)) diss <- d^2/(2*n) else
#     for (i in 1:(n-1))
#       for (j in (i+1):n) 
#         diss[n*(i-1) - i*(i-1)/2 + j-i] <- d[n*(i-1) - i*(i-1)/2 + j-i]^2*wt[i]*wt[j]/(wt[i]+wt[j])
#   return(diss)
# }