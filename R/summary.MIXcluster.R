#' @method summary MIXcluster
#' @title
#'     Summarizing clustering results
#'
#' @description
#'     \code{summary} method for class "\code{MIXcluster}".
#'
#' @param object an object of class "\code{MIXcluster}"
#' @param ...	further arguments passed to or from other methods.
#'
#' @seealso
#' \code{\link{MIXclustering}}
#'
#' @export

summary.MIXcluster <- function(object,...){

  if(class(object)!="MIXcluster") {stop('"object" is not class "MIXcluster"')}
  cat("\n")

  cat("Model based clustering for mixed scale variables.\n")
  cat('\n')

  cat('Call:\n\n')
  cat('\t');print(object$call)
  cat('\n')

  cat('Data description:\n')
  cat('\t Num. rows:\t',object$Y.n,"\n")
  cat('\t Num. variables:\t',object$Y.p,"\n")
  cat('\t Variable types:\n',
      '\t\t',sum(object$Y.var_type=="c"),'continuous\n',
      '\t\t',sum(object$Y.var_type=="o"),'ordinal (categorical, ordered)\n',
      '\t\t',sum(object$Y.var_type=="m"),'nominal (categorical, unordered)\n')
  cat('\n')

  cat('Cluster closest-to-average:\n\n')

  cat('\t Number of clusters:\t',max(object$cluster),"\n\n")

  cat("Cluster:\t rows in data\n")
  for( i in 1:min(10,max(object$cluster)) ) {
    if(sum(object$cluster==i)<=15) {
      cat('\t',i,":\t",paste(which(object$cluster==i),collapse=","),"\n",sep="")
    } else {
      cat('\t',i,":\t",paste(which(object$cluster==i)[1:15],collapse=","),",...(",sum(object$cluster==i)-15,"+)\n",sep="")
    }
  }
  cat("\t(object$cluster for complete information)\n")
  cat('\n')

  cat('\t Averages of Y variables within each cluster:\n')
  print(object$Y.cluster.summary,row.names=F)
}
