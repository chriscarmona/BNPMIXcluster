#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double log_f_post_Omega_cpp( arma::mat Omega,
                             arma::mat Z,
                             arma::mat mu_Z,
                             arma::colvec sampling_prob,
                             arma::mat Lambda ) {

  //////     Metropolis-Hastings for correlations of sigma_Z given by 'Omega'     //////
  // target distribution: log-posterior distribution of 'sigma_jk' //

  double log_post_Omega;

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

  log_post_Omega = 0;
  log_post_Omega += (-0.5) * (n_q+1) * arma::sum( log( abs( Omega.diag() ) ) );
  log_post_Omega += (-0.5) * (n+2-n_q*(n_q-1)) * det(Omega);
  log_post_Omega += (-0.5) * trace( Omega.i() * ( Lambda.i() * S_mat * Lambda.i() ) );

  return log_post_Omega;
}
