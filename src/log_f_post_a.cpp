#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

double log_f_prior_a_cpp( double a,
                          double alpha,
                          double d_0_a,
                          double d_1_a ) {

  // prior distribution of 'a'

  double log_prior_a;
  if( a==0 ) {
    log_prior_a=log(alpha);
  } else {
    log_prior_a=log(1-alpha) + R::dbeta( a, d_0_a, d_1_a, 1 );
  }

  return log_prior_a;

}

// [[Rcpp::export]]
double log_f_post_a_cpp( double a,
                         double b,
                         double alpha,
                         double d_0_a,
                         double d_1_a,
                         arma::colvec mu_star_n_r) {

  // Metropolis-Hastings for 'a' //
  // target distribution: log-posterior distribution of 'a' //

  double log_post_a=0;

  unsigned int r;

  // auxiliar variables
  unsigned int i;
  double aux_double;
  arma::vec aux_vec;

  if( (a<0) | (a>1) ){
    throw std::range_error("The value for a should be in [0,1)");
  }

  //if( (b+a)<0 ){ browser() }
  if( (b+a)<0 ){
    throw std::range_error("The value for b has to be greater than -a");
  }

  r = mu_star_n_r.n_rows;
  if(r>1) {
    aux_vec = arma::linspace<arma::vec>(1, r-1, r-1);
    log_post_a += arma::sum( log( b + a*aux_vec ) );
  }
  aux_double=0;
  for( i=0; i<r ; i++){
    aux_double += lgamma(mu_star_n_r(i)-a)-lgamma(1-a);
  }
  log_post_a += aux_double;
  log_post_a += log_f_prior_a_cpp( a,
                                   alpha,
                                   d_0_a,
                                   d_1_a);

  return log_post_a;
}
