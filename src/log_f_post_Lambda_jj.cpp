#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_f_post_Lambda_jj_cpp( double sigma_jj,
                                 double d_0_z,
                                 double d_1_z,
                                 arma::mat Z,
                                 arma::mat mu_Z,
                                 arma::mat sigma_Z,
                                 arma::colvec sampling_prob ) {

  //////     Metropolis-Hastings for variances of sigma_Z given by 'Lambda'     //////
  // target distribution: log-posterior distribution of 'sigma_jj' //

  double log_f_post_sigma_j;

  unsigned int n;
  unsigned int n_q;

  // auxiliar variables
  unsigned int i;
  arma::mat aux_mat;

  arma::mat S_mat;

  if ( arma::size(Z)!=arma::size(mu_Z) ) {
    throw std::range_error("Problem simulating sigma: Dimension of Z and mu_Z are not equal");
  }

  n = Z.n_rows;
  n_q = Z.n_cols;

  S_mat = arma::zeros<arma::mat>(n_q,n_q);

  for( i=0; i<n; i++ ) {
    aux_mat = Z.row(i) - mu_Z.row(i);
    aux_mat = aux_mat.t();
    S_mat += ( (1/sampling_prob(i)) * (aux_mat * aux_mat.t()) );
  }

  log_f_post_sigma_j = 0;
  log_f_post_sigma_j += - ( d_0_z + n/2 + 1 ) * log(sigma_jj);
  log_f_post_sigma_j += - d_1_z/sigma_jj;
  log_f_post_sigma_j += - trace( sigma_Z.i() * S_mat )/2;

  return log_f_post_sigma_j;

}
