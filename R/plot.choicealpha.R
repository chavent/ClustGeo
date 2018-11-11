#' @title  Plot to choose the mixing parameter
#' @description Plot two curves of explained
#'  inertia (one for \code{D0} and one for \code{D1}) calculated with 
#'  \code{choicealpha}.
#' @param x an object of class \code{choicealpha}.
#' @param norm if TRUE, the normalized explained inertia are plotted. 
#' Otherwise, the explained inertia are plotted.
#' @param lty a vector of size 2 with the line types of the two curves. 
#' See \link{par}
#' @param pch a vector of size 2 specifying the symbol for the points of the 
#' two curves. See \link{par}
#' @param type a vector of size 2 specifying the type of lines of the two 
#' curves. See \link{par}
#' @param col a vector of size 2 specifying the colors the two curves. 
#' See \link{par}
#' @param xlab the title fot the x axis.
#' @param ylab the title fot the y axis.
#' @param legend a vector of size two the the text for the legend of the two curves.
#' @param cex text size in the legend.
#' @param \dots further arguments passed to or from other methods.
#' 
#' @examples
#' data(estuary)
#' D0 <- dist(estuary$dat)
#' D1 <- as.dist(estuary$D.geo) # the geographic distances between the cities
#' range.alpha <- seq(0,1,0.1)
#' K <- 5
#' cr <- choicealpha(D0,D1,range.alpha,K,graph=FALSE)
#' plot(cr,cex=0.8,norm=FALSE,cex.lab=0.8,ylab="pev",
#'          col=3:4,legend=c("socio-demo","geo"), xlab="mixing parameter")
#' plot(cr,cex=0.8,norm=TRUE,cex.lab=0.8,ylab="pev",
#'          col=5:6,pch=5:6,legend=c("socio-demo","geo"), xlab="mixing parameter")
#'          
#' @seealso \code{\link{choicealpha}}
#' 
#' @references 
#' M. Chavent, V. Kuentz-Simonet, A. Labenne, J. Saracco. ClustGeo: an R package
#' for hierarchical clustering with spatial constraints.
#' Comput Stat (2018) 33: 1799-1822.
#' 
#' @export

plot.choicealpha <- function(x,norm=FALSE,lty=1:2,pch=c(8,16),type=c("b","b"),col=1:2,
                             xlab="alpha",ylab=NULL,
                             legend=NULL,cex=1,...)
{
	if (!inherits(x, "choicealpha")) 
      	stop("use only with \"choicealpha\" objects")
  if (is.null(legend))
      legend <- paste("based on",c("D0","D1"))
	if (norm==FALSE)
	{
	  if (is.null(ylab)) 
	    ylab="Q"
	  listpos <- c("topleft", "bottomleft", "topright", "bottomright")
    pos <- listpos[order(c(1-x$Q[1,1],x$Q[1,2],1-x$Q[length(x$range.alpha),2],
                         x$Q[length(x$range.alpha),1]),decreasing=TRUE)[1]]
    graphics::matplot(x$range.alpha,x$Q,ylim=c(0,1), 
            lty=lty,pch=pch,type=type,col=col,
            xlab=xlab,ylab=ylab,...)
    #cex.lab=0.8,cex.main=0.8,cex.axis=0.8
    graphics::legend(pos,legend=legend, col=col,
           lty=lty,pch=pch,bty="n",cex=cex)
	} else
	{
	  if (is.null(ylab)) 
	    ylab="Qnorm"
    listpos <- c("bottomleft","bottomright")
    pos <- listpos[order(c(x$Qnorm[1,2], x$Qnorm[length(x$range.alpha),1]),
                         decreasing=TRUE)[1]]
    graphics::matplot(x$range.alpha,x$Qnorm,ylim=c(0,1), 
            lty=lty,pch=pch,type=type,col=col,
            xlab=xlab,ylab=ylab,...)
    graphics::legend(pos,legend=legend, col=col,
           lty=lty,pch=pch,bty="n",cex=cex)
    graphics::mtext(side=3,paste("of ", round(x$Q[1,1]*100, digits=0),
                                 "%",sep=""),adj=0,cex=cex,...)
    graphics::mtext(side=3,paste("of ",round(x$Q[length(x$range.alpha),2]*100,
                                             digits=0),"%",sep=""),adj=1,cex=cex,...)
  }
}