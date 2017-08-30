
log_f_post_Omega <- function(Omega,
                             Z, mu_Z, sampling_prob,Lambda ) {

  ###     Metropolis-Hastings for correlations of sigma_Z given by 'Omega'     ###
  # target distribution: log-posterior distribution of 'sigma_jk' #

  if(any(dim(Z)!=dim(mu_Z))) {
    cat('\nError: Problem simulating sigma: Dimension of Z and mu_Z are not equal \n')
    stop('Problem simulating sigma: Dimension of Z and mu_Z are not equal')
  }

  n <- nrow(Z)
  n_q <- ncol(Z)

  S_mat <- matrix(0,nrow=n_q,ncol=n_q)
  for(i in 1:nrow(Z)) {
    aux_s <- matrix(Z[i,]-mu_Z[i,],nrow=n_q,ncol=1)
    S_mat <- S_mat + 1/sampling_prob[i] * aux_s %*% t(aux_s)
  }
  #browser()
  log_posterior_Omega <- 0
  log_posterior_Omega <- log_posterior_Omega - (1/2) * (n_q+1) * sum(log(abs(diag(Omega))))
  log_posterior_Omega <- log_posterior_Omega - (1/2) * (n+2-n_q*(n_q-1)) * det(Omega)
  log_posterior_Omega <- log_posterior_Omega - (1/2) * sum(diag( solve(Omega) %*% (solve(Lambda) %*% S_mat %*% solve(Lambda)) ))
  return( log_posterior_Omega )
}

log_f_post_Omega <- compiler::cmpfun(log_f_post_Omega)
