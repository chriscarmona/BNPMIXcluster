#' @importFrom stats rbinom
#' 

# MH Sampling from 'a' #
sampling_a <- function( n_sim_mh=1, a_ini,
                        b, alpha, d_0_a, d_1_a,
                        mu_star_n_r,
                        max_it_time=10*60,n_burn=0,
                        accept_display=F,verbose=F) {

  if( a_ini<0 | a_ini>1 ){
    cat('\nError: the value for "a_ini" has to be in [0,1) \n')
    stop('the value for "a_ini" has to be in [0,1)') }
  if( (b+a_ini)<0 ){
    cat('\nError: the value for "b" has to be greater than -a_ini \n')
    stop('the value for "b" has to be greater than -a_ini')
  }

  # browser()

  n <- sum(mu_star_n_r)
  r <- length(mu_star_n_r)

  # initializing the chain #
  a_chain <- a_ini
  accept_indic <- as.numeric(NA)

  it_t_0 <- Sys.time()
  it_t_i <- as.numeric(Sys.time()-it_t_0)

  while( (length(a_chain)-1 < (n_sim_mh+n_burn)) & (it_t_i<max_it_time) ) {
    if(verbose) {cat(".")}
    it_t_i <- as.numeric(Sys.time()-it_t_0)

    # generate proposal "a_new"
    if(b>0){
      a_prop <- rbinom( n=1, size=1, prob=1-alpha )
      a_prop <- ifelse( a_prop==0, 0, runif(n=1,min=0,max=1) )
    } else {
      a_prop <- runif(n=1,min=max(0,-b),max=1)
    }


    # posterior probability for the current value of a
    log_f_post_a_curr <- log_f_post_a(a=a_chain[length(a_chain)],
                                      b=b,
                                      alpha=alpha,
                                      d_0_a=d_0_a,d_1_a=d_1_a,
                                      mu_star_n_r=mu_star_n_r)

    # posterior probability for the proposal value of a
    log_f_post_a_prop <- log_f_post_a(a=a_prop,
                                      b=b,
                                      alpha=alpha,
                                      d_0_a=d_0_a,d_1_a=d_1_a,
                                      mu_star_n_r=mu_star_n_r)

    log_r <- log_f_post_a_prop - log_f_post_a_curr
    # exp(log_r)

    #if( runif(1,0,1) > lik_ratio ) {
    if(is.na(log_r)) {browser()}
    if( log_r>=0 || runif(1,0,1) < exp(log_r) ) {
      a_chain <- c( a_chain, a_prop )
      accept_indic <- c(accept_indic,1)
      if(verbose) {cat(length(unique(a_chain))-1)}
    } else {
      a_chain <- c( a_chain, a_chain[length(a_chain)] )
      accept_indic <- c(accept_indic,0)
    }
  }

  if(it_t_i>max_it_time) {
    cat('Error: There is a problem simulating from "a" \n')
    stop('There is a problem simulating from "a"')
  }
  #browser()
  if(accept_display) {
    return( list( a_chain = a_chain[(n_burn+2):length(a_chain)],
                  accept_indic = accept_indic[(n_burn+2):length(accept_indic)] ) )
  } else {
    return( a_chain[(n_burn+2):length(a_chain)] )
  }
}
