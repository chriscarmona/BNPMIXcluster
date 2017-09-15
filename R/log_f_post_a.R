
# target distribution: log-posterior distribution of 'a' #
log_f_post_a <- function(a,
                         b,
                         alpha,d_0_a,d_1_a,
                         mu_star_n_r) {

  if( a<0 | a>1 ){
    cat('\nError: The value for "a" has to be in [0,1)\n')
    stop('The value for "a" has to be in [0,1)')
  }
  #if( (b+a)<0 ){ browser() }
  if( (b+a)<0 ){
    cat('\nError: The value for "b" has to be greater than -a \n')
    stop('The value for "b" has to be greater than -a')
  }

  n <- sum(mu_star_n_r)
  r <- length(mu_star_n_r)

  # prior distribution of 'a'
  log_prior_f_a <- function(a,alpha,d_0_a,d_1_a){
    ifelse( a==0 ,
            log(alpha),
            log(1-alpha) + dbeta( x=a, shape1=d_0_a, shape2=d_1_a, log=T )
    ) }

  log_f_post_a <- sum( log(b+a*(1:(r-1)) ) ) + sum( lgamma(x=mu_star_n_r-a)-lgamma(x=1-a) ) + log_prior_f_a(a=a,alpha=alpha,d_0_a=d_0_a,d_1_a=d_1_a)

  return(log_f_post_a)
}
