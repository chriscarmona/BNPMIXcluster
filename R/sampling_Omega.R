
sampling_Omega <- function(n_sim_mh=1,Omega_ini,i,j,delta=4,
                           Z, mu_Z, Lambda, design_prob,
                           max_it_time=10*60, burn_in=0,
                           accept_display=F, verbose=F) {
  # initializing the chain #
  omega_ij_chain <- Omega_ini[i,j]
  accept_indic <- as.numeric(NA)

  it_t_0 <- Sys.time()
  it_t_i <- as.numeric(Sys.time()-it_t_0)

  while( (length(omega_ij_chain)-1 < (n_sim_mh+burn_in)) & (it_t_i<max_it_time) ) {

    if(verbose) {cat(",")}
    it_t_i <- as.numeric(Sys.time()-it_t_0)

    # Generate proposal "omega_ij_new"

    a <- ( det_omega_mod(Omega_ini,1,i=i,j=j) + det_omega_mod(Omega_ini,-1,i=i,j=j) - 2*det_omega_mod(Omega_ini,0,i=i,j=j) ) / 2
    b <- ( det_omega_mod(Omega_ini,1,i=i,j=j) - det_omega_mod(Omega_ini,-1,i=i,j=j) ) / 2
    c <- det_omega_mod(Omega_ini,0,i=i,j=j)

    if(F) {
      # using a grid to sample from Omega
      x <- sort(runif(n=1e4,min=-1,max=1))
      f_x <- a*x^2+b*x+c
      if(!any(f_x>0)) {
        cat('\nError: There is a problem simulating from "Omega"\n')
        stop('There is a problem simulating from "Omega"')
      }
      omega_ij_prop <- sample(x=x[f_x>0],size=1)
      rm(x,f_x)
    }
    if(T) {
      # using a uniform distribution to sample from Omega
      # polinomial roots
      f_r1 <- (-b+sqrt(b^2-4*a*c))/(2*a)
      f_r2 <- (-b-sqrt(b^2-4*a*c))/(2*a)
      if((f_r1<(-1)&f_r2<(-1))|(f_r1>1&f_r2>1)) {
        cat('\nError: There is a problem simulating from "Omega"\n')
        stop('There is a problem simulating from "Omega"')
      }
      sim_last <- omega_ij_chain[length(omega_ij_chain)]
      feasible_int <- sort(c(-1,f_r1,f_r2,1))[c(2,3)]
      aux_o <- diff(feasible_int)/delta

      interval_prop_given_last <- c( max(feasible_int[1],sim_last-aux_o), min(feasible_int[2],sim_last+aux_o) )

      omega_ij_prop <- runif( n=1, min=interval_prop_given_last[1], max=interval_prop_given_last[2] )

      interval_last_given_prop <- c( max(feasible_int[1],omega_ij_prop-aux_o), min(feasible_int[2],omega_ij_prop+aux_o) )

      rm(f_r1,f_r2,sim_last,feasible_int,aux_o)
    }
    Omega_prop <- Omega_ini
    Omega_prop[i,j] <- Omega_prop[j,i] <- omega_ij_prop

    # posterior probability for the current value of omega
    log_f_post_Omega_curr <- log_f_post_Omega( Omega=Omega_ini,
                                               Z=Z, mu_Z=mu_Z, Lambda=Lambda,
                                               design_prob=design_prob )
    #log_f_post_omega_curr

    # posterior probability for the proposal value of omega
    log_f_post_Omega_prop <- log_f_post_Omega( Omega=Omega_prop,
                                               Z=Z, mu_Z=mu_Z, Lambda=Lambda,
                                               design_prob=design_prob )
    # log_f_post_omega_prop

    #log_r <- log_f_post_Omega_prop - log_f_post_Omega_curr
    log_r <- log_f_post_Omega_prop - log_f_post_Omega_curr  + ( log(1/diff(interval_last_given_prop)) - log(1/diff(interval_prop_given_last)) )

    #if(is.nan(log_r)){log_r<--Inf}

    # exp(log_r)
    #browser()
    #if( runif(1,0,1) > lik_ratio ) {
    if( log_r>=0 || runif(1,0,1) < exp(log_r) ) {
      omega_ij_chain <- c( omega_ij_chain, omega_ij_prop )
      accept_indic <- c(accept_indic,1)
      if(verbose) { cat(length(unique(omega_ij_chain))-1) }
    } else {
      omega_ij_chain <- c( omega_ij_chain, omega_ij_chain[length(omega_ij_chain)] )
      accept_indic <- c(accept_indic,0)
    }
  }

  if(it_t_i>max_it_time) {
    cat('\nError: There is a problem simulating from "Sigma"\n')
    stop('There is a problem simulating from "Sigma"')
  }
  # browser()
  if(accept_display) {
    return( list( omega_ij_chain = omega_ij_chain[(burn_in+2):length(omega_ij_chain)],
                  accept_indic = accept_indic[(burn_in+2):length(accept_indic)] ) )
  } else {
    return( omega_ij_chain[(burn_in+2):length(omega_ij_chain)] )
  }
}
