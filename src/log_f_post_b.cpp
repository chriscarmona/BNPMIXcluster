#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_f_post_b_cpp( double b,
                         double a,
                         double d_0_b,
                         double d_1_b,
                         arma::colvec mu_star_n_r ) {

  //////     Metropolis-Hastings for 'b'     //////
  // target distribution: log-posterior distribution of 'b' //
  double log_prior_f_b;
  double log_f_post_b;
  unsigned int n;
  unsigned int r;

  // auxiliar variables
  arma::vec aux_vec;

  n = arma::sum(mu_star_n_r);
  r = mu_star_n_r.n_rows;

  // log-prior distribution of 'b'
  log_prior_f_b = R::dgamma( b+a, d_0_b, d_1_b, 1 );
  
  // log-posterior distribution of 'b'
  log_f_post_b = 0;
  log_f_post_b += log_prior_f_b;
  log_f_post_b += lgamma(b+1)-lgamma(b+n);
  if(r>1) {
    aux_vec = arma::linspace<arma::vec>(1, r-1, r-1);
    log_f_post_b += arma::sum( log( b+a*aux_vec ) );
  }
  
  return log_f_post_b;

}
