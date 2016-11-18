
log_f_post_Lambda_jj <- function(sigma_jj,
                         d_0_z,d_1_z,
                         Z, mu_Z, sigma_Z, sampling_prob
                         ) {

  ###     Metropolis-Hastings for variances of sigma_Z given by 'Lambda'     ###
  # target distribution: log-posterior distribution of 'sigma_jj' #

  if(all(dim(Z)==!dim(mu_Z))) {
    cat('\nError: Problem simulating sigma: Dimension of Z and mu_Z are not equal \n')
    stop('Problem simulating sigma: Dimension of Z and mu_Z are not equal')
  }

  n <- nrow(Z)
  n_q <- ncol(Z)

  S_mat <- matrix(0,nrow=n_q,ncol=n_q)
  for(i in 1:nrow(Z)) {
    aux_s <- matrix(Z[i,]-mu_Z[i,],nrow=n_q,ncol=1)
    S_mat <- S_mat + ( 1/sampling_prob[i] * aux_s %*% t(aux_s) )
  }

  log_f_post_sigma_j <- -(d_0_z+n/2+1)*log(sigma_jj) - d_1_z/sigma_jj - sum(diag( solve(sigma_Z)%*%S_mat ))/2

  return(log_f_post_sigma_j)
}
