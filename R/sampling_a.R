#' Simulation of "\eqn{a}" in the \emph{mixdpcluster} model for bayesian clustering.
#'
#' @description
#' Generates simulation of the posterior distribution of \eqn{a} in the \emph{mixdpcluster} model for bayesian clustering.
#' The simulation is done via Metropolis-Hastings method.
#'
#' @param n number of simulations to generate
#' @param a.ini initialization value
#' @param b parameter \eqn{b} in the posterior distribution of a
#' @param alpha parameter \eqn{\alpha} in the posterior distribution of a
#' @param d_0_a parameter \eqn{d_0^a} in the posterior distribution of a
#' @param d_1_a parameter \eqn{d_1^a} in the posterior distribution of a
#' @param mu_star_n_r bla
#' @param n.burn number of iterations in the simulation considered in the burn-in period.
#' @param n.thin number of iterations discarded between two simulated values (for thinning of the MCMC chain).
#' @param max.time maximum allowed time for the simulation process. The function returns \code{Error} if exceeded.
#' @param verbose if \code{T}, the function reports extra information on progress.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{$a.chain}{A numeric vector with the simulated values from the posterior distribution of \emph{a}}
#'   \item{$accept.indic}{A numeric vector with the simulated values from the posterior distribution of \emph{a}}
#' }
#'
#' @references
#' Carmona C., Nieto-Barajas L., Canale A. (2017). \emph{Model based approach for household clustering with mixed scale variables.}
#'
#' @importFrom stats rbinom
#'
#' @keywords internal

sampling_a <- function( n=1, a.ini,
                        b, alpha, d_0_a, d_1_a,
                        mu_star_n_r,
                        n.burn=0,n.thin=0,
                        max.time=Inf,
                        verbose=F) {

  # Metropolis-Hastings for 'a' #
  # MH Sampling from 'a' #

  if( a.ini<0 | a.ini>1 ){
    cat('\nError: the value for "a.ini" has to be in [0,1) \n')
    stop('the value for "a.ini" has to be in [0,1)') }
  if( (b+a.ini)<0 ){
    cat('\nError: the value for "b" has to be greater than -a.ini \n')
    stop('the value for "b" has to be greater than -a.ini')
  }

  # browser()

  #n.obs <- sum(mu_star_n_r)
  #r <- length(mu_star_n_r)

  # initializing the chain #
  a.chain <- as.numeric(NULL)
  accept.indic <- as.logical(NULL)
  a <- a.ini

  it_t_0 <- Sys.time()
  it_t_i <- as.numeric(Sys.time()-it_t_0)

  # total number of iterations
  n.iter <- n.burn + 1 + (n-1)*(n.thin+1)

  while( ( length(a.chain) < n.iter ) & (it_t_i<max.time) ) {

    if(verbose) {cat(".")}
    it_t_i <- as.numeric(Sys.time()-it_t_0)

    # generate proposal "a_new"
    if(b>0){
      a.prop <- rbinom( n=1, size=1, prob=1-alpha )
      a.prop <- ifelse( a.prop==0, 0, runif(n=1,min=0,max=1) )
    } else {
      a.prop <- runif(n=1,min=max(0,-b),max=1)
    }


    # posterior probability for the current value of a
    log_f_post_a <- log_f_post_a(a=a,
                                      b=b,
                                      alpha=alpha,
                                      d_0_a=d_0_a,d_1_a=d_1_a,
                                      mu_star_n_r=mu_star_n_r)

    # posterior probability for the proposal value of a
    log_f_post_a.prop <- log_f_post_a(a=a.prop,
                                      b=b,
                                      alpha=alpha,
                                      d_0_a=d_0_a,d_1_a=d_1_a,
                                      mu_star_n_r=mu_star_n_r)

    log_r <- log_f_post_a.prop - log_f_post_a
    # exp(log_r)

    #if( runif(1,0,1) > lik_ratio ) {
    if(is.na(log_r)) {browser()}
    if( log_r>=0 || runif(1,0,1) < exp(log_r) ) {
      a <- a.prop
      a.chain <- c( a.chain, a.prop )
      accept.indic <- c(accept.indic,T)
    } else {
      a.chain <- c( a.chain, a )
      accept.indic <- c(accept.indic,F)
    }
  }

  if(it_t_i>max.time) {
    cat('Error: TIMEOUT, There is a problem simulating from "a" \n')
    stop('TIMEOUT, There is a problem simulating from "a"')
  }

  index.iter.out <- seq(from=n.burn+1,to=n.iter,by=(n.thin+1))

  a.sim <- list( a.chain = a.chain[index.iter.out],
                 accept.indic = accept.indic )

  return( a.sim )
}

sampling_a <- compiler::cmpfun(sampling_a)
