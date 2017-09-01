// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

#include <cmath>
#include <R.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// This C++ code implements an Accept/Reject sampler for a single
// Truncated Normal random variable with a mixture of algorithms
// depending on distributional parameters.

/// Check if simpler subalgorithm is appropriate.
inline bool CheckSimple(const double low, ///< lower bound of distribution
                        const double high ///< upper bound of distribution
) {
  // Init Values Used in Inequality of Interest
  double val1 = (2 * sqrt(exp(1.0))) / (low + sqrt(pow(low, 2) + 4));
  double val2 = exp((pow(low, 2) - low * sqrt(pow(low, 2) + 4)) / (4)) ;
  //
  
  // Test if Simple is Preferred
  if (high > low + val1 * val2) {
    return true ;
  } else {
    return false ;
  }
}

/// Draw using algorithm 1.

///
/// Naive Accept-Reject algorithm.
///
inline double UseAlg1(const double low, ///< lower bound of distribution
                      const double high ///< upper bound of distribution
) {
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Init Draw Storage
  double z = 0.0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    z = Rf_rnorm(0.0, 1.0) ;
    
    if (z <= high && z >= low) {
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}

/// Draw using algorithm 2.

///
///  Accept-Reject Algorithm
///

inline double UseAlg2(const double low ///< lower bound of distribution
) {
  // Init Values
  const double alphastar = (low +
                            sqrt(pow(low, 2) + 4.0)
  ) / (2.0) ;
  const double alpha = alphastar ;
  double e = 0 ;
  double z = 0 ;
  double rho = 0 ;
  double u = 0 ;
  //
  
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    e = Rf_rexp(1.0) ;
    z = low + e / alpha ;
    
    rho = exp(-pow(alpha - z, 2) / 2) ;
    u = Rf_runif(0, 1) ;
    if (u <= rho) {
      // Keep Successes
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}

/// Draw using algorithm 3.

///
/// Accept-Reject Algorithm
///

inline double UseAlg3(const double low, ///< lower bound of distribution
                      const double high ///< upper bound of distribution
) {
  // Init Valid Flag
  int valid = 0 ;
  //
  
  // Declare Qtys
  double rho = 0 ;
  double z = 0 ;
  double u = 0 ;
  //
  
  // Loop Until Valid Draw
  while (valid == 0) {
    z = Rf_runif(low, high) ;
    if (0 < low) {
      rho = exp((pow(low, 2) - pow(z, 2)) / 2) ;
    } else if (high < 0) {
      rho = exp((pow(high, 2) - pow(z, 2)) / 2) ;
    } else if (0 < high && low < 0) {
      rho = exp(- pow(z, 2) / 2) ;
    }
    
    u = Rf_runif(0, 1) ;
    if (u <= rho) {
      valid = 1 ;
    }
  }
  //
  
  // Returns
  return z ;
  //
}


/// Draw from an arbitrary truncated normal distribution.

///
/// See Robert (1995): <br />
/// Reference Type: Journal Article <br />
/// Author: Robert, Christian P. <br />
/// Primary Title: Simulation of truncated normal variables <br />
/// Journal Name: Statistics and Computing <br />
/// Cover Date: 1995-06-01 <br />
/// Publisher: Springer Netherlands <br />
/// Issn: 0960-3174 <br />
/// Subject: Mathematics and Statistics <br />
// Start Page: 121 <br />
// End Page: 125 <br />
/// Volume: 5 <br />
/// Issue: 2 <br />
/// Url: http://dx.doi.org/10.1007/BF00143942 <br />
/// Doi: 10.1007/BF00143942 <br />
///

// [[Rcpp::export]]
double rtn1(const double mean,
            const double sd,
            const double low,
            const double high
) {
  // Namespace
  using namespace Rcpp ;
  //
  
  // Init Useful Values
  double draw = 0;
  int type = 0 ;
  int valid = 0 ; // used only when switching to a simplified version
  // of Alg 2 within Type 4 instead of the less
  // efficient Alg 3
  //
  
  // Set Current Distributional Parameters
  const double c_mean = mean ;
  double c_sd = sd ;
  const double c_low = low ;
  const double c_high = high ;
  double c_stdlow = (c_low - c_mean) / c_sd ;
  double c_stdhigh = (c_high - c_mean) / c_sd ; // bounds are standardized
  //
  
  // Map Conceptual Cases to Algorithm Cases
  // Case 1 (Simple Deterministic AR)
  // mu \in [low, high]
  if (0 <= c_stdhigh &&
      0 >= c_stdlow
  ) {
    type = 1 ;
  }
  
  // Case 2 (Robert 2009 AR)
  // mu < low, high = Inf
  if (0 < c_stdlow &&
      c_stdhigh == INFINITY
  ) {
    type = 2 ;
  }
  
  // Case 3 (Robert 2009 AR)
  // high < mu, low = -Inf
  if (0 > c_stdhigh &&
      c_stdlow == -INFINITY
  ) {
    type = 3 ;
  }
  
  // Case 4 (Robert 2009 AR)
  // mu -\in [low, high] & (abs(low) =\= Inf =\= high)
  if ((0 > c_stdhigh || 0 < c_stdlow) &&
      !(c_stdhigh == INFINITY || c_stdlow == -INFINITY)
  ) {
    type = 4 ;
  }
  
  ////////////
  // Type 1 //
  ////////////
  if (type == 1) {
    draw = UseAlg1(c_stdlow, c_stdhigh) ;
  }
  
  ////////////
  // Type 3 //
  ////////////
  if (type == 3) {
    c_stdlow = -1 * c_stdhigh ;
    c_stdhigh = INFINITY ;
    c_sd = -1 * c_sd ; // hack to get two negative signs to cancel out
    
    // Use Algorithm #2 Post-Adjustments
    type = 2 ;
  }
  
  ////////////
  // Type 2 //
  ////////////
  if (type == 2) {
    draw = UseAlg2(c_stdlow) ;
  }
  
  ////////////
  // Type 4 //
  ////////////
  if (type == 4) {
    if (CheckSimple(c_stdlow, c_stdhigh)) {
      while (valid == 0) {
        draw = UseAlg2(c_stdlow) ;
        // use the simple
        // algorithm if it is more
        // efficient
        if (draw <= c_stdhigh) {
          valid = 1 ;
        }
      }
    } else {
      draw = UseAlg3(c_stdlow, c_stdhigh) ; // use the complex
      // algorithm if the simple
      // is less efficient
    }
  }
  
  
  
  // Returns
  return  c_mean + c_sd * draw ;
  //
}


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
        Z(i,n_c+j) = rtn1( mu_Z_ij,
          sqrt(sigma_Z_ij),
          thres_o( cat_Yij ),
          thres_o( cat_Yij+1 ) );
        
        if( ! ( (thres_o( cat_Yij ) < Z(i,n_c+j)) & (Z(i,n_c+j) < thres_o( cat_Yij+1 )) ) ){
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
            Z(i,j_z_s) = rtn1( mu_Z_ij,
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
          Z(i,j_z_s) = rtn1( mu_Z_ij,
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
              Z(i,j_z_s) = rtn1( mu_Z_ij,
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
