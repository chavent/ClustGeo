#' @title Empirical choice of the mixing parameter
#' @description This function calculates the proportion (resp. normalized proportion) of explained inertia of
#'  the partitions in K clusters obtained with the Ward-like \code{hclustgeo} procedure for a range 
#'  of mixing parameters alpha. When the proportion (resp. normalized proportion) of explained inertia based on D0 decreases,
#'   the proportion (resp. normalized proportion)
#'  of explained inertia based on D1 increases. The plot of these criteria can help the user in the choice of the
#'  mixing parameter alpha.
#' @param D0 an object of class "dist" with the dissimilarities between the n observations. 
#' The function \code{\link{as.dist}} can be used to transform an object of class matrix to object of class "dist".
#' @param D1 an object of class "dist" with other dissimilarities between the same n observations. 
#' @param range.alpha  a vector of real values between 0 and 1. 
#' @param K  the  number of clusters. 
#' @param wt vector with the weights of the observations. By default, wt=NULL corresponds to
#' the case where all observations are weighted by 1/n.
#' @param scale if TRUE the two dissimilarity matrix are scaled i.e. divided by their max.
#' @param graph if TRUE the two graphics (proportion and normalized proportion of explained inertia) are drawn.
# @keywords internal
#' @export
#' @examples
#' data(estuary)
#' D0 <- dist(estuary$dat) # the socio-demographic distances
#' D1 <- as.dist(estuary$D.geo) # the geographic distances between the cities
#' range.alpha <- seq(0,1,0.1)
#' K <- 5
#' cr <- choicealpha(D0,D1,range.alpha,K,graph=TRUE)
#' cr$Q # proportion of explained pseudo inertia
#' cr$Qnorm # normalized proportion of explained pseudo inertia
#' @references 
#' M.chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco.  ClustGeo:  an R package 
#' for hierarchical clustering with spatial constraints	arXiv:1707.03897 [stat.CO]


choicealpha <- function(D0, D1, range.alpha,K, wt=NULL,scale=TRUE,graph=TRUE) {

  if (is.null(D0)) stop("D0 must be an argument",call.=FALSE)
  if (is.null(D1)) stop("D1 must be an argument",call.=FALSE)
  if (class(D0)!="dist")
    stop("DO must be of class dist (use as.dist)",call.=FALSE)
  if (class(D1)!="dist")
    stop("D1 must be of class dist (use as.dist)",call.=FALSE)
  
  n.alpha<-length(range.alpha)
  n <- as.integer(attr(D1, "Size"))
  if (is.null(n)) 
    stop("invalid dissimilarities",call.=FALSE)
  if (is.na(n) || n > 65536L) 
    stop("size cannot be NA nor exceed 65536",call.=FALSE)
  if (n < 2) 
    stop("must have n >= 2 objects to cluster",call.=FALSE)
  if (!is.null(D1) && length(D0) != length(D1))
    stop("the two dissimilarity structures must have the same size",call.=FALSE)
  if ((max(range.alpha) >1) || (max(range.alpha) < 0))
    stop("Values range.alpha must be in [0,1]",call. = FALSE)
  
  if (scale==TRUE) {
    #scaled dissimilarity matrices 
    D0 <- D0/max(D0)
    D1 <- D1/max(D1)
  }
  
  if (is.null(wt)) 
    wt <- rep(1/n, n)
  
  # within-cluster inertia obtained either with D0 or D1
  W <- matrix(0,length(range.alpha),2)
  rownames(W)  <- paste("alpha=", range.alpha, sep="")
  colnames(W) <- c("W0","W1")
  for (i in 1:length(range.alpha)) {
    tree <- hclustgeo(D0,D1,range.alpha[i],scale=scale,wt=wt)
    part <- cutree(tree,k=K)
    W[i,1] <- withindiss(D0,part,wt)
    W[i,2] <- withindiss(D1,part,wt)
  }
  
  # total inertia obtained with either with D0 or D1
  T0 <-  inertdiss(D0,wt=wt)
  T1 <-  inertdiss(D1,wt=wt) 
   
  # proportion of explained inertia obtained either with D0 orD1 
  Q <- matrix(0,length(range.alpha),2)
  rownames(Q)  <- rownames(W)
  colnames(Q) <- c("Q0","Q1")
  Q[,1] <- 1-W[,1]/T0
  Q[,2] <- 1-W[,2]/T1
  
  # normalized proportion of explained inertia 
  Qnorm <- matrix(0,length(range.alpha),2)
  rownames(Qnorm)  <- rownames(W)
  colnames(Qnorm) <- c("Q0norm","Q1norm")
  Qnorm[,1] <- Q[,1]/Q[1,1]
  Qnorm[,2] <- Q[,2]/Q[length(range.alpha),2]
  
  
  if (graph==TRUE)
  {
    listpos <- c("topleft","bottomleft","topright","bottomright")
    pos <- listpos[order(c(1-Q[1,1],Q[1,2],1-Q[length(range.alpha),2],
                         Q[length(range.alpha),1]),decreasing=TRUE)[1]]
    matplot(range.alpha,Q,xlab="alpha",ylim=c(0,1), 
            ylab="Q",type="b",pch=c(8,16),lty=1:2,
            main=paste("K=",K,"clusters"))
    #cex.lab=0.8,cex.main=0.8,cex.axis=0.8
    legend(pos,legend=paste("based on",c("D0","D1")), col=1:2,lty=1,pch=16,bty="n",cex=1)
    
    listpos <- c("bottomleft","bottomright")
    pos <- listpos[order(c(Qnorm[1,2],Qnorm[length(range.alpha),1]),decreasing=TRUE)[1]]
    matplot(range.alpha,Qnorm,xlab="alpha",ylim=c(0,1), 
            ylab="Qnorm",type="b",pch=c(8,16),lty=1:2,
            main=paste("K=",K,"clusters"))
    legend(pos,legend=paste("based on",c("D0","D1")), col=1:2,lty=1,pch=16,bty="n",cex=1)
    mtext(side=3,paste("of ",round(Q[1,1]*100,digits=0),"%",sep=""),cex=1,adj=0)
    mtext(side=3,paste("of ",round(Q[length(range.alpha),2]*100,digits=0),"%",sep=""),cex=1,adj=1)
  }
  retlist <- list(Q=Q,Qnorm=Qnorm,range.alpha=range.alpha,K=K)
  class(retlist) <- "choicealpha"
  return(retlist)
}
