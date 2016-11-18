#' @method print MIXcluster
#' @export

print.MIXcluster <- function(x,...){

  if(class(x)!="MIXcluster") {stop('x is not class "MIXcluster"')}
  cat("\n")

  cat("Model based clustering for mixed scale variables.\n")

  cat('\n')
  cat('Call:\n\n')
  print(x$call)
  cat('\n')

  cat('Cluster closest-to-average:\n\n')

  cat('\t Number of clusters:\t',max(x$cluster),"\n\n")

  cat("Cluster:\t rows in data\n")
  for( i in 1:min(10,max(x$cluster)) ) {
    if(sum(x$cluster==i)<=15) {
      cat('\t',i,":\t",paste(which(x$cluster==i),collapse=","),"\n",sep="")
    } else {
      cat('\t',i,":\t",paste(which(x$cluster==i)[1:15],collapse=","),",...(",sum(x$cluster==i)-15,"+)\n",sep="")
    }
  }
  cat("\t(x$cluster for complete information)\n")
  cat('\n')

}
