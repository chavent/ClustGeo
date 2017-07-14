#' @title Hierarchical clustering with geographical contraints
#' @description This function implements a Ward-like hierarchical clustering algorithm 
#' including soft contiguity constraints. This algorithm takes as input two dissimilarity matrices 
#' D0 and D1 and a mixing parameter alpha between 0 an 1. The dissimilarities can be non 
#' euclidean and the weights of the observations can be non uniform. 
#' The first matrix gives the dissimilarities in the "feature space" 
#' (socio-demographic variables or grey levels for instance). The second matrix gives the dissimilarities in the 
#'  "constraint" space. For instance, D1 can be a matrix of geographical distances or a matrix 
#' build from the contiguity matrix C.
#' The mixing parameter alpha sets the importance of the constraint in the clustering procedure. 
#' @param D0 an object of class "dist" with the dissimilarities between the n observations. 
#' The function \code{\link{as.dist}} can be used to transform an object of class matrix to object of class "dist".
#' @param D1 an object of class "dist" with other dissimilarities between the same n observations. 
#' @param alpha  a real value between 0 and 1. This mixing parameter gives the relative importance of D0 compared to D1.
#' By default, this parameter is equal to 0 and D0 is used alone in the clustering process.  
#' @param wt vector with the weights of the observations. By default, wt=NULL corresponds to
#' the case where all observations are weighted by 1/n.
#' @param scale if TRUE the two dissimilarity matrix D0 and D1 are scaled i.e. divided by their max. If D1=NULL, 
#' this parameter is no used and D0 is not scaled.
#' @return Returns an object of class \code{\link{hclust}}.
#' @details The criterion minimized at each stage is a convex combination of the homogeneity criterion 
#'  calculated with D0 and the homogeneity criterion calculated with D1. 
#'  The parameter alpha (the weight of this convex combination) controls the weight of the constraint in the quality of the solutions. 
#'  When alpha increases, the homogeneity calculated with D0 decreases whereas the homogeneity calculated with D1 increases. 
#' @examples
#' data(estuary)
#' # with one dissimilarity matrix
#' w <- estuary$map@data$POPULATION # non uniform weights 
#' D <- dist(estuary$dat)
#' tree <- hclustgeo(D,wt=w)
#' sum(tree$height)
#' inertdiss(D,wt=w)
#' inert(estuary$dat,w=w)
#' plot(tree,labels=FALSE)
#' part <- cutree(tree,k=5)
#' sp::plot(estuary$map,border="grey",col=part)
#' 
#' # with two dissimilarity matrix
#' D0 <- dist(estuary$dat) # the socio-demographic distances
#' D1 <- as.dist(estuary$D.geo) # the geographical distances
#' alpha <- 0.2 # the mixing parameter
#' tree <- hclustgeo(D0,D1,alpha=alpha,wt=w)
#' plot(tree,labels=FALSE)
#' part <- cutree(tree,k=5)
#' sp::plot(estuary$map,border="grey",col=part)
#' @references 
#' M.chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco.  ClustGeo:  an R package 
#' for hierarchical clustering with spatial constraints	arXiv:1707.03897 [stat.CO]
#' @export


hclustgeo <- function(D0, D1=NULL, alpha=0,scale=TRUE, wt=NULL) {
  if (class(D0)!="dist")
    stop("DO must be of class dist (use as.dist)",call.=FALSE)
  if (!is.null(D1) && (class(D1)!="dist"))
    stop("D1 must be of class dist (use as.dist)",call.=FALSE)
  n <- as.integer(attr(D0, "Size"))
  if (is.null(n)) 
    stop("invalid dissimilarities",call.=FALSE)
  if (is.na(n) || n > 65536L) 
    stop("size cannot be NA nor exceed 65536",call.=FALSE)
  if (n < 2) 
    stop("must have n >= 2 objects to cluster",call.=FALSE)
  if (!is.null(D1) && length(D0) != length(D1))
    stop("the two dissimilarity structures must have the same size",call.=FALSE)
  if ((max(alpha) >1) || (max(alpha) < 0))
    stop("Values alpha must be in [0,1]",call. = FALSE)
  
  if ((scale==TRUE) && (!is.null(D1))) {
    #Normalized dissimilarity matrices 
    D0 <- D0/max(D0)
    D1 <- D1/max(D1)
  }
  #Dissimilarity measure (aggregation measure) between singletons
  delta0 <- wardinit(D0,wt)
  if (!is.null(D1)) 
    delta1 <- wardinit(D1,wt) else delta1 <- 0
  delta <-(1-alpha)*delta0 + alpha*delta1
  
  #Hierarchical clustering
  res <- hclust(delta,method="ward.D",members=wt)
  return(res)
}


#' estuary data
#' @name estuary 
#' @format  The R dataset estuary is a list of three objects: 
#' \itemize{
#' \item{dat: a data frame with the description of the n=303 municipalities on p=4 socio-demographic variables.}
#'  \item{D.geo: a matrix with the geographical distances between the town hall of the n=303 municipalities.}
#'  \item{map: an object of class \code{SpatialPolygonsDataFrame} with the map of the gironde estuary.}
#'  }
#' @source Original data are issued from the French population census of National Institute 
#' of Statistics and Economic Studies for year 2009. The agricultural surface has been 
#' calculated on data coming from the French National Institute of Geographical and Forestry 
#' Information. The calculation of the ratio and recoding of categories have been made by 
#' Irstea Bordeaux.
#' @description Data refering to n=303 french municipalities  of gironde estuary (a south-ouest french county).
#' The data are issued from the French population census conducted by the National Institute 
#' of Statistics and Economic Studies. The dataset is an extraction of four quantitative 
#' socio-economic variables for a subsample of 303 french municipalities located on the
#' atlantic coast between Royan and Mimizan. \code{employ.rate.city} is the employment rate 
#' of the municipality, that is the ratio of the number of individuals who have a job to 
#' the population of working age (generally defined, for the purposes of international 
#' comparison, as persons of between 15 and 64 years of age). \code{graduate.rate} refers 
#' to the level of education of the population that is the highest degree declared by the 
#' individual. It is defined here as the ratio for the whole population having completed 
#' a diploma equivalent or of upper level to two years of higher education 
#' (DUT, BTS, DEUG, nursing and social training courses, license, maitrise, master, DEA, DESS, doctorate, or Grande Ecole diploma). 
#' \code{housing.appart} is the ratio of apartment housing. \code{agri.land} is the part of 
#' agricultural area of the municipality.
#' @keywords data
#' @references 
#' M.chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco.  ClustGeo:  an R package 
#' for hierarchical clustering with spatial constraints	arXiv:1707.03897 [stat.CO]
#' @examples
#' data(estuary)
#' names(estuary)
#' head(estuary$dat)
#' sp::plot(estuary$map)
NULL

#' @importFrom graphics legend matplot mtext
NULL 
#' @importFrom stats as.dist cutree hclust
NULL
#' @importMethodsFrom sp plot
NULL
#' @importMethodsFrom sp coordinates
NULL
#' @importClassesFrom sp SpatialPolygonsDataFrame
NULL
#' @importFrom spdep nb2mat
NULL
#' @importFrom spdep poly2nb
NULL

