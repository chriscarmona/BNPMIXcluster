#' @method plot MIXcluster
#'
#' @title
#'     Plotting clustering results for "\code{MIXcluster}" objects
#'
#' @description
#'     Plotting method for objects inheriting from class "\code{MIXcluster}".
#'
#' @param x an object of class "\code{MIXcluster}"
#' @param type what type of plot should be drawn. Possible types are:
#'        \describe{
#'            \item{\code{"heatmap"}}{ (default) draws a heatmap of the average similarity matrix for the effective iterations of the MCMC.}
#'            \item{\code{"chain"}}{ for the evolution and histograms of the chains for parameters in the model.  }
#'            }
#' @param chain.obj if \code{type="chain"}, this specifies what chain will be plotted. Possible types are:
#'        \describe{
#'            \item{\code{"n.cluster"}}{ (default) for the number of clusters.}
#'            \item{\code{"a"}}{ for the \eqn{a} parameter of the model.}
#'            \item{\code{"b"}}{ for the \eqn{b} parameter of the model.}
#'            \item{\code{"Lambda"}}{ one plot for each element in the diagonal of the \eqn{\Lambda} matrix of the model (variance of latent variables).}
#'            \item{\code{"Omega"}}{ one plot for each element above the diagonal of the \eqn{\Omega} matrix of the model (correlation between latent variables).}
#'            \item{\code{"all"}}{ for all of the above.}
#'            }
#' @param ...	further arguments passed to or from other methods.
#'
#' @seealso
#' \code{\link{MIXclustering}}
#'
#' @importFrom graphics abline barplot hist layout par plot
#' @importFrom stats heatmap
#'
#' @export

plot.MIXcluster <- function(x,
                            type=c("heatmap","chain")[1],
                            chain.obj=c("n.cluster","a","b","Lambda","Omega","all")[1],
                            ...
                            ) {

  if(is.element("heatmap",type)) {
    heatmap(1-x$cluster.matrix.avg)
  }

  if( is.element("chain",type) ) {

    plot_mcmc <- function(x,...){
      #op <- par()
      hist_n <- hist(x,plot=F)
      layout(matrix(1:2,nrow=1,ncol=2),widths=c(4/5,1/5))
      par(mar=c(4,4,2,0.25))
      plot( x,type="l",lwd=2,col="darkblue",...); abline(h=0,col="grey70")
      par(mar=c(4,0.25,2,2))
      barplot( hist_n$counts,axes=F,space=0,xlim=c(0,max(hist_n$counts)),horiz=T)
      #par(op)
    }

    if( (is.element("all",chain.obj)|is.element("n.cluster",chain.obj)) ) {
      plot_mcmc(x$MC.values$n.clusters,main="Number of clusters",ylab="Number of clusters",xlab="iteration")
    }
    if( (is.element("all",chain.obj)|is.element("a",chain.obj)) ) {
      plot_mcmc(x$MC.values$a,main="a",ylab="a",xlab="iteration")
    }
    if( (is.element("all",chain.obj)|is.element("b",chain.obj)) ) {
      plot_mcmc(x$MC.values$b,main="b",ylab="b",xlab="iteration")
    }
    if( (is.element("all",chain.obj)|is.element("Omega",chain.obj)) ) {
      for(i in 1:(dim(x$MC.values$Omega)[1]-1)){
        for(j in (i+1):dim(x$MC.values$Omega)[2]){
          plot_mcmc(x$MC.values$Omega[i,j,],main=paste("Omega[",i,",",j,"]",sep=""),ylab=paste("Omega[",i,",",j,"]",sep=""),xlab="iteration")
        }
      }
      rm(i,j)
    }
    if( (is.element("all",chain.obj)|is.element("Lambda",chain.obj)) ) {
      for(j in 1:ncol(x$MC.values$Lambda)){
        plot_mcmc(x$MC.values$Lambda[,j],main=paste("Lambda[",j,",",j,"]",sep=""),ylab=paste("Lambda[",j,",",j,"]",sep=""),xlab="iteration")
      }
      rm(j)
    }
  }
}
