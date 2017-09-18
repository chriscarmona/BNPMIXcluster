
sampling_b <- function( n_sim_mh=1, b_ini,
                        a, d_0_b, d_1_b, eta=1,
                        mu_star_n_r,
                        max.time=10*60, n.burn=0,
                        accept_display=T, verbose=F,
                        USING_CPP=TRUE ) {

  ###     Metropolis-Hastings for 'b'     ###
  # MH Sampling from 'b' #

  n <- sum(mu_star_n_r)
  r <- length(mu_star_n_r)

  if( a<0 | a>1 ){
    cat('\nError: The value for "a" has to be in [0,1) \n')
    stop('The value for "a" has to be in [0,1)')
  }
  if( (b_ini+a)<0 ){
    cat('\nError: The initial value for "b" has to be greater than -a \n')
    stop('Error: The initial value for "b" has to be greater than -a') }

  # initializing the chain #
  b_chain <- b_ini
  accept_indic <- as.numeric(NA)

  it_t_0 <- Sys.time()
  it_t_i <- as.numeric(Sys.time()-it_t_0)

  while( (length(b_chain)-1 < (n_sim_mh+n.burn)) & (it_t_i<max.time) ) {
    if(verbose) {cat(".")}
    it_t_i <- as.numeric(Sys.time()-it_t_0)

    # generate proposal "b_new"

    # wrong!
    #b_prop <- b_chain[length(b_chain)] + rnorm(1,0,eta)
    #while( (b_prop-a)<0 ) { b_prop <- b_chain[length(b_chain)] + rnorm(1,0,eta) }

    # correct
    b_prop_interval <- c(max(-a,b_chain[length(b_chain)]-eta),b_chain[length(b_chain)]+eta)
    b_prop <- runif(1,b_prop_interval[1],b_prop_interval[2])

    # posterior probability for the current and proposed value of b
    if( USING_CPP ) {
      log_f_post_b_curr <- log_f_post_b_cpp(b=b_chain[length(b_chain)],
                                            a=a,
                                            d_0_b=d_0_b,d_1_b=d_1_b,
                                            mu_star_n_r=mu_star_n_r)
      log_f_post_b_prop <- log_f_post_b_cpp(b=b_prop,
                                            a=a,
                                            d_0_b=d_0_b,d_1_b=d_1_b,
                                            mu_star_n_r=mu_star_n_r)
    } else {
      log_f_post_b_curr <- log_f_post_b( b = b_chain[length(b_chain)],
                                         a = a,
                                         d_0_b = d_0_b,
                                         d_1_b = d_1_b,
                                         mu_star_n_r = mu_star_n_r )
      log_f_post_b_prop <- log_f_post_b( b = b_prop,
                                         a = a,
                                         d_0_b = d_0_b,
                                         d_1_b = d_1_b,
                                         mu_star_n_r = mu_star_n_r )
    }
    
    b_cur_interval_given_prop <- c(max(-a,b_prop-eta),b_prop+eta)

    log_r <- (log_f_post_b_prop - log_f_post_b_curr) + ( log(1/diff(b_cur_interval_given_prop)) - log(1/diff(b_prop_interval)) )
    # exp(log_r)

    #if( runif(1,0,1) > lik_ratio ) {
    if( log_r>=0 || runif(1,0,1) < exp(log_r) ) {
      b_chain <- c( b_chain, b_prop )
      accept_indic <- c(accept_indic,1)
      if(verbose) {cat(length(unique(b_chain))-1)}
    } else {
      b_chain <- c( b_chain, b_chain[length(b_chain)] )
      accept_indic <- c(accept_indic,0)
    }
  }

  if(it_t_i>max.time) {
    cat('\nError: There is a problem simulating from "sigma_jk" \n')
    stop('There is a problem simulating from "sigma_jk"')
  }
  #browser()
  if(accept_display) {
    return( list( b_chain = b_chain[(n.burn+2):length(b_chain)],
                  accept_indic = accept_indic[(n.burn+2):length(accept_indic)] ) )
  } else {
    return( b_chain[(n.burn+2):length(b_chain)] )
  }
}

sampling_b <- compiler::cmpfun(sampling_b)
