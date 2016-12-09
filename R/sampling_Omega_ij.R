#' Simulation of "\eqn{\Omega_{i,j}}" in the \emph{mixdpcluster} model for bayesian clustering.
#'
#' @description
#' Generates simulation of the posterior distribution of the \eqn{(i,j)} element of the \eqn{\Omega} matrix in the \emph{mixdpcluster} model for bayesian clustering.
#' The simulation is done via Metropolis-Hastings method.
#'
#' @param n number of simulations to be generated
#' @param Omega.ini matrix \eqn{\Omega} with an initialization value for \eqn{\Omega_{i,j}}.
#' @param i indicates the row for \eqn{\Omega_{i,j}}
#' @param j indicates the column for \eqn{\Omega_{i,j}}
#' @param delta defines the maximum jump on each iteration of the MCMC as \eqn{1/delta} of the feasible interval for \eqn{\Omega_{i,j}}
#' @param Z bla
#' @param mu_Z bla
#' @param Lambda bla
#' @param sampling_prob bla
#' @param n.burn number of iterations in the simulation considered in the burn-in period.
#' @param n.thin number of iterations discarded between two simulated values (for thinning of the MCMC chain).
#' @param max.time maximum allowed time for the simulation process. The function returns \code{Error} if exceeded.
#' @param verbose if \code{T}, the function reports extra information on progress.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{$omega_ij.chain}{A numeric vector with the simulated values from the posterior distribution of \eqn{\Omega_{i,j} }. }
#'   \item{$accept.indic}{A logical vector indicating whether or not \code{omega_ij.prop} was accepted.}
#' }
#'
#' @references
#' Carmona C., Nieto-Barajas L., Canale A. (2017). \emph{Model based approach for household clustering with mixed scale variables.}
#'
#' @keywords internal

sampling_Omega_ij <- function(n=1,Omega.ini,i,j,delta=4,
                           Z, mu_Z, Lambda, sampling_prob,
                           n.burn=0,n.thin=0,
                           max.time=Inf,
                           verbose=F) {

  ###     Metropolis-Hastings for correlations of sigma_Z given by 'Omega'     ###
  # MH Sampling from 'omega' #

  det_omega_mod <- function(X,x_ij.new,i,j) {
    X[i,j] <- X[j,i] <- x_ij.new
    return(det(X))
  }

  # measuring execution time #
  it_t_0 <- Sys.time()
  it_t_i <- as.numeric(Sys.time()-it_t_0)

  # feasible interval for omega_ij
  # such that Omega is positive definite
  a <- ( det_omega_mod(Omega.ini,1,i=i,j=j) + det_omega_mod(Omega.ini,-1,i=i,j=j) - 2*det_omega_mod(Omega.ini,0,i=i,j=j) ) / 2
  b <- ( det_omega_mod(Omega.ini,1,i=i,j=j) - det_omega_mod(Omega.ini,-1,i=i,j=j) ) / 2
  c <- det_omega_mod(Omega.ini,0,i=i,j=j)
  f_r1 <- (-b+sqrt(b^2-4*a*c))/(2*a)
  f_r2 <- (-b-sqrt(b^2-4*a*c))/(2*a)
  feasible_int <- sort(c(-1,f_r1,f_r2,1))[c(2,3)]
  if((f_r1<(-1)&f_r2<(-1))|(f_r1>1&f_r2>1)) {
    cat('\nError: There is a problem simulating from "Omega"\n')
    stop('There is a problem simulating from "Omega"')
  }


  # initializing the chain #
  omega_ij.chain <- as.numeric(NULL)
  accept.indic <- as.logical(NULL)
  Omega <- Omega.ini

  # total number of iterations
  n.iter <- n.burn + 1 + (n-1)*(n.thin+1)

  while( ( length(omega_ij.chain) < n.iter ) & (it_t_i<max.time) ) {

    if(verbose) {cat(".")}
    it_t_i <- as.numeric(Sys.time()-it_t_0)

    # Generate proposal "omega_ij_new"
    # using a uniform distribution to sample from Omega
    # polinomial roots


    max.jump <- diff(feasible_int)/delta

    interval_prop_given_last <- c( max(feasible_int[1],Omega[i,j]-max.jump), min(feasible_int[2],Omega[i,j]+max.jump) )

    omega_ij.prop <- runif( n=1, min=interval_prop_given_last[1], max=interval_prop_given_last[2] )

    interval_last_given_prop <- c( max(feasible_int[1],omega_ij.prop-max.jump), min(feasible_int[2],omega_ij.prop+max.jump) )

    # Omega matrix with omega_ij.prop
    Omega.prop <- Omega.ini
    Omega.prop[i,j] <- Omega.prop[j,i] <- omega_ij.prop

    # posterior probability for the current value of omega
    log_f_post_Omega <- log_f_post_Omega( Omega=Omega,
                                          Z=Z, mu_Z=mu_Z, Lambda=Lambda,
                                          sampling_prob=sampling_prob )

    # posterior probability for the proposal value of omega
    log_f_post_Omega.prop <- log_f_post_Omega( Omega=Omega.prop,
                                               Z=Z, mu_Z=mu_Z, Lambda=Lambda,
                                               sampling_prob=sampling_prob )

    log_r <- log_f_post_Omega.prop - log_f_post_Omega  + ( log(1/diff(interval_last_given_prop)) - log(1/diff(interval_prop_given_last)) )


    if( log_r>=0 || runif(1,0,1) < exp(log_r) ) {
      Omega <- Omega.prop
      omega_ij.chain <- c( omega_ij.chain, omega_ij.prop )
      accept.indic <- c(accept.indic,T)
    } else {
      omega_ij.chain <- c( omega_ij.chain, Omega[i,j] )
      accept.indic <- c(accept.indic,F)
    }
  }

  if(it_t_i>max.time) {
    cat('\nError: TIMEOUT, There is a problem simulating from "Omega"\n')
    stop('TIMEOUT, There is a problem simulating from "Omega"')
  }

  index.iter.out <- seq(from=n.burn+1,to=n.iter,by=(n.thin+1))

  omega_ij_sim <- list( omega_ij.chain=omega_ij.chain[index.iter.out],
                        accept.indic=accept.indic )

  return( omega_ij_sim )

}

sampling_Omega_ij <- compiler::cmpfun(sampling_Omega_ij)
