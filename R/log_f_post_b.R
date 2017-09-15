
# target distribution: log-posterior distribution of 'b' #
log_f_post_b <- function(b,
                         a,
                         d_0_b,d_1_b,
                         mu_star_n_r) {

  n <- sum(mu_star_n_r)
  r <- length(mu_star_n_r)

  # log-prior distribution of 'b' given 'a'
  log_prior_f_b <- dgamma( b+a, shape=d_0_b, rate=d_1_b, log=T )
  log_f_post_b <- lgamma(x=b+1)-lgamma(x=b+n) + sum( log( b+a*(1:(r-1)) ) ) + log_prior_f_b

  return(log_f_post_b)
}
