#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
#include "../inst/include/BNPMIXcluster.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat get_latents_cpp( arma::mat Y,
                           arma::ucolvec var_type,
                           arma::mat mu_Z,
                           arma::mat sigma_Z,
                           arma::mat Z_ini,
                           bool verbose=false ) {

  ////// This function simulates the latent variables for a specified Y //////
  //require(truncnorm)

  ///// INPUT /////
  // Y : asummed to be in order: continuous, then ordinal, then nominal

  ///// OUTPUT /////
  // Z : matrix with values of latent variables

  // number of observations in Y
  unsigned int n_obs;
  n_obs = Y.n_rows;

  // number of variables in Y
  unsigned int p;
  p = Y.n_cols;

  // auxiliar variables
  unsigned int i=0;
  unsigned int j=0;
  unsigned int s=0;
  arma::ucolvec aux_ucolvec;
  arma::vec aux_vec;
  arma::colvec cat_Yj;
  unsigned int cat_Yij;

  // possible variable classes that are allowed
  arma::vec var_type_all;
  var_type_all << 1 << 2 << 3 << arma::endr;

  // checking input consistency
  if( var_type.n_rows != p ) {
    throw std::range_error("The number of columns in Y have to be equal to the lenght of vector var_type");
  }

  if( arma::any(var_type<1) || arma::any(var_type>3) ) {
    throw std::range_error("Elements in var_type are not supported");
  }

  // number of variables by type
  unsigned int n_c;
  unsigned int n_o;
  unsigned int n_m;
  n_c = arma::sum( var_type==var_type_all(0) );
  n_o = arma::sum( var_type==var_type_all(1) );
  n_m = arma::sum( var_type==var_type_all(2) );

  if(verbose) {
    Rcpp::Rcout << "n_c=" << n_c << ", n_o=" << n_o << ", n_m=" << n_m << std::endl;
  }

  // TRUE
  // p==n_c+n_o+n_m;

  // vector with the number of latent variables that will be used for each categorical variable
  arma::ucolvec n_m_l;
  arma::ucolvec n_m_l_cum;

  if( n_m>0 ) {
    n_m_l = arma::zeros<arma::ucolvec>(n_m);

    // which columns of Y are nominal?
    aux_ucolvec = find( var_type==var_type_all(2) );

    for ( j=0; j<n_m; j++ ) {

      aux_vec = Y.col(n_c+n_o+j);
      aux_ucolvec = arma::find_unique( aux_vec, false);
      cat_Yj = aux_vec(aux_ucolvec);

      // number of categories in that ordinal variable
      // k-1 latents are needed for a categorical variable with k categories
      n_m_l(j) = cat_Yj.n_rows - 1;
    }
  } else {
    n_m_l = arma::zeros<arma::ucolvec>(1);
  }

  // Compute the total number of latents that will be needed //
  unsigned int n_q;
  n_m_l_cum = arma::cumsum(n_m_l);
  n_m_l_cum.insert_rows(0,1);
  n_q = n_c + n_o + arma::sum(n_m_l);

  if(verbose) {
    Rcpp::Rcout << "n_m_l= ";
    for(i=0;i<n_m_l.n_rows;i++) {
      Rcpp::Rcout << n_m_l(i) << " ";
    }
    Rcpp::Rcout << std::endl;
  }

  // Matriz Z with latent variables
  arma::mat Z = arma::zeros<arma::mat>( n_obs, n_q );

  if(verbose) {
    Rcpp::Rcout << "n_obs=" << n_obs << ", n_q=" << n_q << std::endl;
  }

  ////////// Simulating Latent values //////////

  //////   Latent: Continuous   //////
  // Continuous variables will remain the same values //
  if( n_c > 0 ) {
    if(verbose) {
      Rcpp::Rcout << "Generating continuous latents..." << std::endl;
    }
    // which columns of Y are continuous?
    aux_ucolvec = find( var_type==var_type_all(0) );

    if( size(Z_ini)==size(arma::zeros<arma::mat>(1,1)) ) {
      // standardize the continuos variables
      for( j=0; j<n_c; j++ ) {
        Z.col(j) = ( Y.col(j) - arma::mean( Y.col(j) ) )/arma::stddev( Y.col(j) );
      }
    } else {
      for( j=0; j<n_c; j++ ) {
        Z.col(j) = Z_ini.col(j);
      }
    }
  }

  //////   Latent: Categorical   //////

  double mu_Z_ij;
  double sigma_Z_ij;
  arma::mat Sigma_12;
  arma::mat Sigma_22;
  arma::mat Z_i_not_j;
  arma::mat mu_Z_i_not_j;

  // number of categories in each categorical variable
  double K_o;

  ////   Latent: Categorical, Ordinal   ////
  if( n_o>0 ) {

    if(verbose) {
      Rcpp::Rcout << "Generating ordinal latents..." << std::endl;
    }

    // vector of thresholds for each ordinal variable
    arma::vec thres_o;

    for( j=0; j<n_o; j++ ) {

      // categories in that ordinal variable
      aux_vec = Y.col(n_c+j);
      aux_ucolvec = arma::find_unique( aux_vec, false);
      cat_Yj = aux_vec(aux_ucolvec);

      // number of categories in that ordinal variable
      K_o = cat_Yj.n_rows;

      // obtain the thresholds dividing (-Inf,Inf) into K_o[i] intervals
      // thresholds of length 4 with a fixed variance of 1
      thres_o = arma::linspace<arma::vec>(-floor((K_o-2)/2), ceil((K_o-2)/2), K_o-1);
      thres_o *= 4;
      thres_o.insert_rows(0,1);
      thres_o.insert_rows(thres_o.n_rows,1);
      thres_o(0)=-INFINITY;
      thres_o(thres_o.n_rows-1)=INFINITY;

      // simulating ordinal latents //

      // j goes one for each ordinal in Y,  which is the same as each ordinal latent in Z!

      for( i=0; i<n_obs; i++ ) {
        aux_ucolvec = arma::find( cat_Yj==Y(i,n_c+j), 1 );
        cat_Yij = aux_ucolvec(0,0);

        if( size(Z_ini)==size(arma::zeros<arma::mat>(1,1)) ) {
          // if there is no given initial value for Z[i,n_c+j]
          // the simulated ordinal latent Z[i,n_c+j]
          // come from a distribution with mean mu_Z[i,n_c+j] and variance sigma_Z[n_c+j,n_c+j]
          mu_Z_ij = mu_Z(i,n_c+j);
          sigma_Z_ij = sigma_Z(n_c+j,n_c+j);
        } else {
          // if an initial value for Z was given
          Sigma_12 = sigma_Z.row(n_c+j); Sigma_12.shed_col(n_c+j);

          Sigma_22 = sigma_Z;
          Sigma_22.shed_row(n_c+j); Sigma_22.shed_col(n_c+j);

          Z_i_not_j = Z_ini.row(i); Z_i_not_j.shed_col(n_c+j);
          mu_Z_i_not_j = mu_Z.row(i); mu_Z_i_not_j.shed_col(n_c+j);

          aux_vec = Sigma_12 * Sigma_22.i() * ( Z_i_not_j.t() - mu_Z_i_not_j.t() );
          mu_Z_ij = mu_Z( i, n_c+j ) + aux_vec(0,0);

          aux_vec = Sigma_12 * Sigma_22.i() * Sigma_12.t();
          sigma_Z_ij = sigma_Z( n_c+j , n_c+j ) + aux_vec(0,0);
        }

        // Sampling from truncated normal
        Z(i,n_c+j) = BNPMIXcluster::rtn1( mu_Z_ij,
          sqrt(sigma_Z_ij),
          thres_o( cat_Yij ),
          thres_o( cat_Yij+1 ) );

        if( ! ( thres_o( cat_Yij ) < Z(i,n_c+j) & Z(i,n_c+j) < thres_o( cat_Yij+1 ) ) ){
          throw std::range_error("There was a problem simulating an ordinal latent Z_ij");
        }
      }
    }

  }

  //////   Latent: Categorical, Nominal   //////

  if( n_m>0 ) {
    if(verbose) {
      Rcpp::Rcout << "Generating nominal latents..." << std::endl;
    }
    unsigned int j_z_s;
    unsigned int j_z_max;

    for( j=0; j<n_m; j++ ) {
      // categories in that nominal variable
      aux_vec = Y.col(n_c+n_o+j);
      aux_ucolvec = arma::find_unique( aux_vec, false);
      cat_Yj = aux_vec(aux_ucolvec);

      for ( i=0 ; i<n_obs; i++ ) {
        // simulating one by one
        aux_ucolvec = arma::find( cat_Yj==Y(i,n_c+n_o+j), 1 );
        cat_Yij = aux_ucolvec(0,0);

        if( cat_Yij == cat_Yj( cat_Yj.n_rows-1 ) ) {
          // observation Y[i,j] IS in the last class

          // in this case:
          // all values of latents have to be negative

          for( s=0; s<cat_Yj.n_rows-1; s++ ) {
            j_z_s = n_c + n_o + n_m_l_cum(j) + s; // column of z corresponding with this value

            if( size(Z_ini)==size(arma::zeros<arma::mat>(1,1)) ) {
              // if there is no given initial value for Z(i,j_z_s)
              // the simulated latent Z(i,j_z_s)
              // come from a distribution with mean mu_Z(i,j_z_s) and variance sigma_Z( j_z_s , j_z_s )
              mu_Z_ij = mu_Z( i , j_z_s );
              sigma_Z_ij = sigma_Z( j_z_s , j_z_s );
            } else {
              // if an initial value for Z was given
              Sigma_12 = sigma_Z.row(j_z_s); Sigma_12.shed_col(j_z_s);

              Sigma_22 = sigma_Z;
              Sigma_22.shed_row(j_z_s); Sigma_22.shed_col(j_z_s);

              Z_i_not_j = Z_ini.row(i); Z_i_not_j.shed_col(j_z_s);
              mu_Z_i_not_j = mu_Z.row(i); mu_Z_i_not_j.shed_col(j_z_s);

              aux_vec = Sigma_12 * Sigma_22.i() * ( Z_i_not_j.t() - mu_Z_i_not_j.t() );
              mu_Z_ij = mu_Z( i, j_z_s ) + aux_vec(0,0);

              aux_vec = Sigma_12 * Sigma_22.i() * Sigma_12.t();
              sigma_Z_ij = sigma_Z( j_z_s , j_z_s ) + aux_vec(0,0);
            }
            // simulates negative gaussian numbers //
            Z(i,j_z_s) = BNPMIXcluster::rtn1( mu_Z_ij,
              sqrt(sigma_Z_ij),
              -INFINITY,
              0 );
          }

        } else {
          // observation Y[i,j] IS NOT in the last class

          // in this case:
          // at least one latent >0
          // and the maximum value of the associated Z's corresponds to the value of the class in Y

          ////// first, simulates the maximum //////

          // column of Z[i,] corresponding with the category of Y[i,]
          // it will take the maximum value of the corresponding latent
          j_z_max = n_c + n_o + n_m_l_cum(j) + cat_Yij;
          j_z_s = j_z_max;

          if( size(Z_ini)==size(arma::zeros<arma::mat>(1,1)) ) {
            // if there is no given initial value for Z(i,j_z_s)
            // the simulated latent Z(i,j_z_s)
            // come from a distribution with mean mu_Z(i,j_z_s) and variance sigma_Z( j_z_s , j_z_s )
            mu_Z_ij = mu_Z( i , j_z_s );
            sigma_Z_ij = sigma_Z( j_z_s , j_z_s );
          } else {
            // if an initial value for Z was given
            Sigma_12 = sigma_Z.row(j_z_s); Sigma_12.shed_col(j_z_s);

            Sigma_22 = sigma_Z;
            Sigma_22.shed_row(j_z_s); Sigma_22.shed_col(j_z_s);

            Z_i_not_j = Z_ini.row(i); Z_i_not_j.shed_col(j_z_s);
            mu_Z_i_not_j = mu_Z.row(i); mu_Z_i_not_j.shed_col(j_z_s);

            aux_vec = Sigma_12 * Sigma_22.i() * ( Z_i_not_j.t() - mu_Z_i_not_j.t() );
            mu_Z_ij = mu_Z( i, j_z_s ) + aux_vec(0,0);

            aux_vec = Sigma_12 * Sigma_22.i() * Sigma_12.t();
            sigma_Z_ij = sigma_Z( j_z_s , j_z_s ) + aux_vec(0,0);
          }

          // simulates gaussian number //
          Z(i,j_z_s) = BNPMIXcluster::rtn1( mu_Z_ij,
            sqrt(sigma_Z_ij),
            0,
            INFINITY );

          // then, simulates the rest //
          for( s=0; s<cat_Yj.n_rows-1; s++ ) {
            j_z_s = n_c + n_o + n_m_l_cum(j) + s; // column of z corresponding with this value
            if (j_z_s!=j_z_max) {
              if( size(Z_ini)==size(arma::zeros<arma::mat>(1,1)) ) {
                // if there is no given initial value for Z(i,j_z_s)
                // the simulated latent Z(i,j_z_s)
                // come from a distribution with mean mu_Z(i,j_z_s) and variance sigma_Z( j_z_s , j_z_s )
                mu_Z_ij = mu_Z( i , j_z_s );
                sigma_Z_ij = sigma_Z( j_z_s , j_z_s );
              } else {
                // if an initial value for Z was given
                Sigma_12 = sigma_Z.row(j_z_s); Sigma_12.shed_col(j_z_s);

                Sigma_22 = sigma_Z;
                Sigma_22.shed_row(j_z_s); Sigma_22.shed_col(j_z_s);

                Z_i_not_j = Z_ini.row(i); Z_i_not_j.shed_col(j_z_s);
                mu_Z_i_not_j = mu_Z.row(i); mu_Z_i_not_j.shed_col(j_z_s);

                aux_vec = Sigma_12 * Sigma_22.i() * ( Z_i_not_j.t() - mu_Z_i_not_j.t() );
                mu_Z_ij = mu_Z( i, j_z_s ) + aux_vec(0,0);

                aux_vec = Sigma_12 * Sigma_22.i() * Sigma_12.t();
                sigma_Z_ij = sigma_Z( j_z_s , j_z_s ) + aux_vec(0,0);
              }

              // simulates gaussian number //
              Z(i,j_z_s) = BNPMIXcluster::rtn1( mu_Z_ij,
                sqrt(sigma_Z_ij),
                -INFINITY,
                Z(i,j_z_max) );
            }
          }
        }
      }
    }

  }
  return Z;
}
