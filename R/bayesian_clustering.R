
bayesian_clustering <- function( Y,
                                 data_class=c("c","o","m")[rep(1,ncol(data))],
                                 expansion_f=NULL,
                                 design_prob=NULL,
                                 n_iter=2000,
                                 burn_in=0,
                                 thin_k=1,
                                 a_fix=NULL,b_fix=NULL,
                                 alpha=0,
                                 d_0_z=1, d_1_z=1, kappa=5,
                                 d_0_a=1, d_1_a=1,
                                 d_0_b=1, d_1_b=1, eta=2, delta=4,
                                 d_0_mu=1, d_1_mu=1,
                                 Z_sim_file=NULL,
                                 get_param_evol=F,
                                 max_it_time=30*60,
                                 cluster_file=NULL,
                                 review_MH_chain=F
                                 ) {

  #browser()

  #require(MASS) # ginv(sigma_Z)
  #require(truncnorm)

  dev_verbose<-F

  # Y <- data.frame( as.matrix(Y) )

  # possible variable classes that are allowed
  var_classes <- c("c","o","m")

  # checking input consistency
  if( length(data_class) != ncol(Y) ) {
    cat('\nError: The number of columns in "Y" have to be equal to the lenght of vector "data_class"\n')
    stop('The number of columns in "Y" have to be equal to the lenght of vector "data_class"')
  }
  if( any(!is.element(data_class,var_classes)) ) {
    cat('\nError: Elements in "data_class" have to be one of ',paste(var_classes,collapse = ","),'\n')
    stop('Elements in "data_class" have to be one of ',paste(var_classes,collapse = ","))
  }

  # Sorting Y columns
  Y_new_order <- c( which(is.element(data_class,var_classes[1])),
                    which(is.element(data_class,var_classes[2])),
                    which(is.element(data_class,var_classes[3])))
  data_class <- data_class[Y_new_order]
  Y <- Y[,Y_new_order]
  rm(Y_new_order)

  n <- nrow(Y)
  p <- ncol(Y)

  # number of variables by type
  n_c <- sum( is.element( data_class, var_classes[1] ) )
  n_o <- sum( is.element( data_class, var_classes[2] ) )
  n_m <- sum( is.element( data_class, var_classes[3] ) )

  # changes the colnames of Y for simplicity and standarization
  colnames(Y) <- paste("var_",
                       c(rep(var_classes[1],n_c),rep(var_classes[2],n_o),rep(var_classes[3],n_m)),"_",
                       formatC( unlist(mapply(seq,1,c(n_c,n_o,n_m),length.out=c(n_c,n_o,n_m))) , width = 2, flag = '0'),
                       sep="")

  if(is.null(design_prob)) {

    if(is.null(expansion_f)) {
      expansion_f<-rep(1,n)
    } else {
      if(length(expansion_f)!=n) {
        cat('\nError: There is an inconsistency between "expansion_f" and the number of rows in "Y"\n')
        stop('There is an inconsistency between "expansion_f" and the number of rows in "Y"')
      }
    }

    design_prob <- 1/expansion_f # Sample design probabilities



    if(length(design_prob)!=n) {
      cat('\nError: There is an inconsistency between "design_prob" and the number of rows in "Y"\n')
      stop('There is an inconsistency between "design_prob" and the number of rows in "Y"')
    }
  }

  #design_prob <- design_prob/sum(design_prob) # should add up 1
  #browser()
  if(n_c>0){
    ### scale continuos variables to have mean 0 and sd 1 ###
    Y[,1:n_c] <- scale(Y[,1:n_c],center=T,scale=T)
  }

  #####     Simulating latent variables 'Z' from 'Y'     #####

  latents_info <- get_latents( Y, data_class )

  Z <- latents_info$Z

  #browser()

  n_c <- latents_info$n_c
  n_o <- latents_info$n_o
  n_m_l <- latents_info$n_m_l
  n_q <- n_c+n_o+sum(n_m_l)

  #n_q==latents_info$n_q # TRUE

  K_o <- latents_info$K_o

  rm(latents_info); gc()

  if( !is.null(Z_sim_file) ) {
    write.csv(Y,file="Y.csv",quote=F,na='',row.names=F)
    write.csv(Z,file=gsub(".csv",paste("_00.csv",sep=""),Z_sim_file),quote=F,na="",row.names=F)
  }



  ##### Missing data #####
  # Only consider rows with complete information, i.e. no missing data allowed
  aux_na_Y <- apply(!is.na(Z),1,all)
  Z <- Z[aux_na_Y,]
  if(sum(!aux_na_Y)>0){
    cat(sum(!aux_na_Y),' rows with missing data were removed\n')
  }
  # vector for returning Z to Y incorporating back the NAs
  map_z_to_Y <-cumsum(aux_na_Y)
  map_z_to_Y[duplicated(map_z_to_Y)]<-NA
  #rm(aux_na_Y)

  n <- dim(Z)[1]

  ##### Start: Initializing the chain #####

  # sigma_Z: variances in Lambda
  Lambda <- diag(apply(Z,2,sd))

  # unit variances for "sigma_Z"
  #browser()
  aux_var1_Z <- rep(F,n_q)
  if( n_o>0 ){
    # unit variance for binary ordinal variables
    # aux_var1_Z[n_c+which(K_o==2)] <- T
    # unit variance for ALL ordinal variables
    aux_var1_Z[(n_c+1):(n_c+n_o)] <- T
  }
  if( n_m>0 ){
    # unit variance for categorical (not ordinal) variables
    aux_var1_Z[(n_c+n_o+1):(n_c+n_o+sum(n_m_l))] <- T
  }
  diag(Lambda)[aux_var1_Z] <- 1
  #rm(aux_var1_Z)

  #cat("aux_var1_Z: ",paste(aux_var1_Z,collapse=", "),"\n")

  # sigma_Z: correlations in Omega
  Omega <- cor(Z)

  # sigma_Z
  sigma_Z <- Lambda %*% Omega %*% Lambda

  # checking that sigma_Z[j,k] <= sigma_Z[j,j] * sigma_Z[k,k]
  aux_cond_sigma <- sigma_Z < outer(diag(sigma_Z),diag(sigma_Z))
  diag(aux_cond_sigma) <- T
  if ( !all(aux_cond_sigma) ) {
    cat('\nError: It must happen that: sigma_Z[j,k] < sigma_Z[j,j] * sigma_Z[k,k] \n')
    stop('It must happen that: sigma_Z[j,k] < sigma_Z[j,j] * sigma_Z[k,k] ')
  }
  rm(aux_cond_sigma)

  # mu_Z <- Z
  # it is not necessary to define an object for mu_Z, since it is is fully characterized by
  # its unique values 'mu_star'
  # and its mapping 'mu_star_map'

  # matrix with unique values of mu_Z
  mu_star <- Z[!duplicated(Z[,1:2]),]
  # mapping each mu_Z to mu_star
  mu_star_map <- match(data.frame(t(Z[,1:2])), data.frame(t(mu_star[,1:2]))) # only matches with the two first values for computational efficiency
  # how many repeated values of each mu_star
  mu_star_n_r <- as.numeric(table(mu_star_map))

  #rm(mu_Z); gc()

  # 'sigma_mu'
  sigma_mu <- diag( diag(sigma_Z) )

  if( (nrow(mu_star)!=length(mu_star_n_r)) | (nrow(Z)!=length(mu_star_map)) ) {
    cat('Error: there is an inconsistency between "mu_star", "mu_star_n_r" and "mu_star_map"\n')
    stop('there is an inconsistency between "mu_star", "mu_star_n_r" and "mu_star_map"')
  }

  a <- ifelse( is.null(a_fix), 0, a_fix )

  b <- ifelse( is.null(b_fix), 0.1, b_fix )

  ##### End: Initializing the chain #####

  #####
  #mu_Z_sim <- array(NA,dim=c(n,n_q,n_iter))
  #sigma_Z_sim <- array(NA,dim=c(n_q,n_q,n_iter))
  mu_star_map_sim <- data.frame(iter_0=1:n)
  #colnames(mu_star_map_sim) <- paste("iter_",1:ncol(mu_star_map_sim),sep="")
  #rownames(mu_star_map_sim) <- NULL

  #mu_star_n_r_sim <- list()
  #mu_star_probs <- array(dim=c(n,n+1,n_iter))

  ##### Monitoring acceptance rate for MH #####

  if(review_MH_chain) {

    Lambda_sim <- matrix(as.numeric(NA),nrow=n_iter,ncol=ncol(sigma_Z) )
    Omega_sim <- array(as.numeric(NA),dim=c(nrow(sigma_Z),ncol(sigma_Z),n_iter))
    a_sim <- as.numeric(NULL)
    b_sim <- as.numeric(NULL)

    Lambda_accept <- matrix(as.numeric(NA),nrow=n_iter,ncol=ncol(sigma_Z) )
    Omega_accept <- array(as.numeric(NA),dim=c(nrow(sigma_Z),ncol(sigma_Z),n_iter))
    a_accept <- as.numeric(NULL)
    b_accept <- as.numeric(NULL)

  } else {

    Lambda_sim <- NULL
    Omega_sim <- NULL
    a_sim <- NULL
    b_sim <- NULL

    Lambda_accept <- NULL
    Omega_accept <- NULL
    a_accept <- NULL
    b_accept <- NULL

  }

  #if(dev_verbose) {
  cat("*** Clustering estimation started ***\n")
  #cat('     Sys.time() = ',as.character(Sys.time()),sep='')
  #}
  time_sim <- Sys.time()

  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)

  lsa_iter <- c(ls(),"iter_i","lsa_iter")

  for( iter_i in 1:n_iter) {
    # iter_i<-1
    clock_iter_i <- Sys.time()
    if(dev_verbose) {
      cat('\n**-- iter_i=',iter_i,' --**\n\n')
    }


    ### (a) Sampling mu_Z ###

    #mu_Z_new <- mu_Z ; mu_Z_new[] <- NA
    if(dev_verbose) {
      cat('Sampling "mu_Z": ',nrow(Z),' values...\n',sep='')
    }
    # Creating a copy to keep control of original values of mu_star
    mu_star_new <- mu_star
    mu_star_map_new <- mu_star_map
    mu_star_n_r_new <- mu_star_n_r
    for( i in 1:nrow(Z) ) {
      # i<-1
      if(dev_verbose) {
        if (i%%10==0) { cat( i , ", " ) }
      }
      # getting the probabilities #
      v_i <- solve( solve(design_prob[i] * sigma_Z) + solve(sigma_mu) )
      # u_i <- v_i %*% solve(design_prob[i] * sigma_Z) %*% matrix(Z[i,],nrow=n_q,ncol=1)
      u_i <- v_i %*% solve(design_prob[i] * sigma_Z) %*% matrix(mu_star[mu_star_map[i],],nrow=n_q,ncol=1)

      # eliminates the non-simmetry due to numerical precision
      if( max( abs(v_i-t(v_i)) )<1e-7 ) {v_i <- v_i - (v_i-t(v_i))/2} else {browser()}

      # number of diferent values of mu_star[-1]
      mu_star_n_r_temp <- mu_star_n_r
      mu_star_n_r_temp[mu_star_map[i]] <- mu_star_n_r_temp[mu_star_map[i]]-1
      r_i <- sum(mu_star_n_r_temp>0)

      D_0 <- ( b + a * r_i ) * mvtnorm::dmvnorm( x=Z[i,] , mean=rep(0,n_q), sigma=design_prob[i]*sigma_Z+sigma_mu )
      D_values <- as.numeric(NULL)
      for(r in 1:length(mu_star_n_r_temp)) {
        D_values[r] <- ( mu_star_n_r_temp[r] - a ) * mvtnorm::dmvnorm( x=Z[i,] , mean=mu_star[r,], sigma=design_prob[i]*sigma_Z )
      }
      D_values[D_values<0] <- 0

      if(sum( D_0, D_values )==0) {
        #browser()
        cat('\nError: There is a problem with "D_r" in the simulation of mu_Z \n: sum( D_0, D_values )==0\n')
        stop('There is a problem with "D_r" in the simulation of mu_Z \n: sum( D_0, D_values )==0')
      }

      # Calculating probabilities
      prob_0 <- D_0 / sum( D_0, D_values )
      prob_values <- D_values / sum( D_0, D_values )

      #mu_star_probs[i,1+0:length(prob_values),iter_i] <- c(prob_0,prob_values)

      # simulating mu #
      mu_class <- which( is.element( rmultinom(n=1,size=1,prob=c(prob_0,prob_values)), 1 ) )
      if(mu_class==1) {
        # generates a new value for mu_Z
        #browser()
        #if(iter_i==5) {browser()}
        mu_star_new <- rbind( mu_star_new, mvtnorm::rmvnorm( n=1 , mean=u_i, sigma=v_i ) )
        mu_star_map_new[i] <- nrow(mu_star_new)
        mu_star_n_r_new[mu_star_map[i]] <- mu_star_n_r_new[mu_star_map[i]]-1
        mu_star_n_r_new <- c(mu_star_n_r_new,1)

      } else {
        # get an existing value
        #cat(i,mu_class-1,"\n")
        mu_star_n_r_new[mu_star_map[i]] <- mu_star_n_r_new[mu_star_map[i]]-1
        mu_star_n_r_new[mu_class-1] <- mu_star_n_r_new[mu_class-1]+1
        mu_star_map_new[i] <- mu_class-1
      }
      if(any(is.na(mu_star_new))){browser()}
    } #finished simulating from mu_Z

    # Identifies values of mu_star that have no Z elements pointing to them
    aux_no_map <- is.element(mu_star_n_r_new,0)

    # Replaces new simulated values of mu_star
    mu_star <- mu_star_new[!aux_no_map,]
    mu_star_map <- mu_star_map_new - cumsum(aux_no_map)[mu_star_map_new]
    mu_star_n_r <- mu_star_n_r_new[!aux_no_map]

    rm(mu_star_new,mu_star_map_new,mu_star_n_r_new,mu_star_n_r_temp,
       v_i,u_i,D_0,D_values,prob_0,prob_values,aux_no_map,
       r,r_i,mu_class,
       i); gc()
    if(dev_verbose) {
      cat('...Done! \n')
    }

    ### (b) Sampling mu_star ###

    # Creating a copy to keep control of original values of mu_star
    mu_star_new <- mu_star
    mu_star_new[] <- NA
    if(dev_verbose) {
      cat('Sampling "mu_star": ',nrow(mu_star),' values...\n',sep='')
    }
    sigma_mu_inv <- solve(sigma_mu)
    sigma_Z_inv <- solve(sigma_Z)
    for( j in 1:nrow(mu_star) ) {
      # j<-1
      if(dev_verbose) {
        if (j%%10==0) { cat( j , ", " ) }
      }
      I_j <- is.element(mu_star_map,j)

      if(T){
        v_i_star <- solve(sum(1/design_prob[I_j]) * sigma_Z_inv + sigma_mu_inv )
      } else {
        for( k in 1:sum(I_j) ) {
          if(k==1) {
            v_i_star <- (1/design_prob[I_j])[k] * sigma_Z_inv + sigma_mu_inv
          } else {
            v_i_star <- v_i_star + (1/design_prob[I_j])[k] * sigma_Z_inv + sigma_mu_inv
          }
        }
        v_i_star <- solve(v_i_star)
        rm(k)
      }
      #browser()
      weighted_z <- apply( matrix(1/design_prob[I_j],nrow=sum(I_j),ncol=n_q,byrow=F) * matrix(Z[I_j,],nrow=sum(I_j),ncol=n_q), 2, sum )
      weighted_z <- matrix(weighted_z,nrow=n_q,ncol=1)
      u_i_star <- v_i_star %*% solve(sigma_Z) %*% weighted_z

      mu_star_new[j,] <- mvtnorm::rmvnorm(1,mean=u_i_star,sigma=v_i_star)
    }
    if(dev_verbose) {
      cat('...Done! \n')
    }
    # Replaces new simulated values of mu_star
    mu_star <- mu_star_new

    rm(mu_star_new,
       v_i_star,weighted_z,u_i_star,
       I_j,sigma_mu_inv,sigma_Z_inv,
       j); gc()


    ### (c) Sampling "sigma_mu" ###
    if(dev_verbose) {
      cat( 'Sampling "sigma_mu":\n' )
    }
    sigma_mu_new <- sigma_mu
    sigma_mu_new[] <- 0; diag(sigma_mu_new) <- NA
    for( j in 1:n_q ) {
      if(dev_verbose) {
        cat(j,", ")
      }
      shape_gamma <- d_0_mu+length(mu_star_n_r)/2
      rate_gamma <- d_1_mu+(1/2)*sum(mu_star[,j]^2)
      sigma_mu_new[j,j] <- 1/rgamma( n=1, shape=shape_gamma, rate=rate_gamma )
    }
    sigma_mu <- sigma_mu_new
    rm(shape_gamma,rate_gamma); gc()
    if(dev_verbose) {
      cat('...Done! \n')
    }

    ### Sampling sigma_Z ###

    if( !matrixcalc::is.positive.definite(sigma_Z) ) {
      cat("*****\nProcess finished because 'sigma_Z' is not positive definite!\n*****");
      return()
    }

    sigma_Z_new <- sigma_Z
    Lambda_new <- Lambda
    Omega_new <- Omega

    #rm(sigma_Z,Lambda,Omega)

    ### (d) Sampling sigma_Z, variances in Lambda ###
    if(dev_verbose) {
      cat( 'Sampling sigma_Z: variances in Lambda \n' )
    }


    for( j_sigma in which(!aux_var1_Z)) { # all elements of sigma_Z with variance different of 1
    #for( j_sigma in 1:n_q) { # all elements of sigma_Z with variance different of 1

      if(dev_verbose) {
        cat(j_sigma,", ")
      }

      aux_Lambda <-  sampling_sigma_jj( n_sim_mh=1, sigma_jj_ini=Lambda_new[j_sigma,j_sigma]^2,j=j_sigma,
                                        d_0_z=d_0_z, d_1_z=d_1_z, kappa=kappa,
                                        Z=Z, mu_Z=mu_star[mu_star_map,], sigma_Z=sigma_Z_new, design_prob=design_prob,
                                        max_it_time=10*60,burn_in=0,
                                        accept_display=review_MH_chain, verbose=F)

      if(review_MH_chain) {
        Lambda_new[j_sigma,j_sigma] <- sqrt(aux_Lambda[[1]])
        Lambda_accept[iter_i,j_sigma] <- aux_Lambda[[2]]
      } else {
        Lambda_new[j_sigma,j_sigma] <- sqrt(aux_Lambda[[1]])
        #hist(aux_Lambda,50)
        #abline(v=Lambda_new[j_sigma,j_sigma]^2,col="red")
      }

      #browser()

      # Element with unitary variance in Lambda_new
      # diag(Lambda_new)[aux_var1_Z] <- 1

      # updates 'sigma_Z' after updating each element of 'Lambda'
      # because 'Lambda' sampling IS dependent of 'sigma_Z'
      sigma_Z_new <- Lambda_new %*% Omega_new %*% Lambda_new

      # eliminates the non-simmetry due to numerical precision
      if( max( abs(sigma_Z_new-t(sigma_Z_new)) ) < 1e-7 ) {sigma_Z_new <- sigma_Z_new - (sigma_Z_new-t(sigma_Z_new))/2} else {browser()}

      if( !matrixcalc::is.positive.definite(sigma_Z_new) ) {
        cat("*****\nProcess finished because 'sigma_Z_new' is not positive definite!\n*****");
        return()
      }
      rm(aux_Lambda)
    }

    rm(j_sigma); gc()

    if(dev_verbose) {
      cat('...Done! \n')
    }


    ### (e) Sampling sigma_Z, correlations in Omega ###
    if(dev_verbose) {
      cat( 'Sampling sigma_Z: correlations in Omega \n' )
    }

    if( !matrixcalc::is.positive.definite(Omega_new) ) {
      cat("*****\nProcess finished because 'Omega_new' is not positive definite!\n*****");
      return()
    }

    for(i_omega in 2:dim(Omega_new)[1]) {
      for(j_omega in 1:(i_omega-1) ) {
        if(dev_verbose) {
          cat('(',i_omega,',',j_omega,') ',sep='')
        }
        aux_omega_ij_new <- sampling_Omega(n_sim_mh=1,Omega_ini=Omega_new,i=i_omega,j=j_omega,
                                           delta=delta,
                                           Z=Z, mu_Z=mu_star[mu_star_map,], Lambda=Lambda_new, design_prob=design_prob,
                                           max_it_time=10*60, burn_in=0,
                                           accept_display=review_MH_chain, verbose=F)

        if(review_MH_chain) {
          omega_ij_new <- aux_omega_ij_new[[1]]
          Omega_accept[i_omega,j_omega,iter_i] <- Omega_accept[j_omega,i_omega,iter_i] <- aux_omega_ij_new[[2]]
        } else {
          #hist(aux_omega_ij_new,50)
          #abline(v=Omega_new[i_omega,j_omega],col="red")
          omega_ij_new <- aux_omega_ij_new
        }

        Omega_new[i_omega,j_omega] <- Omega_new[j_omega,i_omega] <- omega_ij_new
        if(dev_verbose) {
          if(j_omega==(i_omega-1)) {cat('\n')}
        }
        if( !matrixcalc::is.positive.definite(Omega_new) ) {
          cat("     Process finished because 'Omega_new' is not positive definite!\n");
          return()
        }

      }
    }
    rm(i_omega,j_omega,aux_omega_ij_new,omega_ij_new);gc()

    # updates 'sigma_Z' after all element in 'Omega' are sampled
    # only once because 'Omega' sampling IS NOT dependent of 'sigma_Z'
    # browser()
    sigma_Z_new <- Lambda_new %*% Omega_new %*% Lambda_new

    if( max( abs(sigma_Z_new-t(sigma_Z_new)) ) < 1e-7 ) {sigma_Z_new <- sigma_Z_new - (sigma_Z_new-t(sigma_Z_new))/2} else {browser()}

    if( !matrixcalc::is.positive.definite(sigma_Z_new) ) {
      cat("*****\nProcess finished because 'sigma_Z_new' is not positive definite!\n*****");
      return()
    }

    #browser()

    if(F) {
      # Sampling sigma_Z in the old "tricky" way

      # sigma_Z: variances in Lambda
      if(dev_verbose) {
        cat('Sampling sigma_Z[j,j]...\n')
      }
      for( j in 1:n_q) {
        if(dev_verbose) {
          cat(j,", ")
        }
        # j<-1
        shape_gamma <- d_0_z + (n/2)
        rate_gamma <- d_1_z + (1/2) * sum( (1/design_prob) * ( Z[,j] - mu_star[mu_star_map,j] )^2 )
        Lambda_new[j,j] <- 1/rgamma(n=1,shape=shape_gamma,rate=rate_gamma)
      }
      rm(shape_gamma,rate_gamma); gc()

      diag(Lambda_new)[aux_var1_Z] <- 1
      #rm(aux_var1_Z)

      # sigma_Z: correlations in Omega
      Omega_new <- cor(Z)

      sigma_Z_new <- Lambda_new %*% Omega_new %*% Lambda_new

    }

    # eliminates the non-simmetry due to numerical precision
    if( max( abs(sigma_Z_new-t(sigma_Z_new)) ) < 1e-7 ) {sigma_Z_new <- sigma_Z_new - (sigma_Z_new-t(sigma_Z_new))/2} else {browser()}

    #if( !matrixcalc::is.positive.definite(sigma_Z) ) {browser()}
    if( !matrixcalc::is.positive.definite(sigma_Z_new) ) {
      cat("*****\nProcess finished because 'sigma_Z' is not positive definite!\n*****");
      return()
    }

    sigma_Z <- sigma_Z_new
    Lambda <- Lambda_new
    Omega <- Omega_new
    if(dev_verbose) {
      cat('...Done! \n')
    }

    # (f) Sampling "a"
    if(dev_verbose) {
      cat( 'Sampling "a":\n' )
    }
    if(is.null(a_fix)){
      aux_a_new <- sampling_a( n_sim_mh=1, a_ini=a,
                           b=b, alpha=alpha, d_0_a=d_0_a, d_1_a=d_1_a,
                           mu_star_n_r=mu_star_n_r,
                           max_it_time=max_it_time, burn_in=0,
                           accept_display=review_MH_chain, verbose=F )

      #plot(a_new)
      #hist(a_new)

      if(review_MH_chain) {
        a_new <- aux_a_new[[1]]
        a_accept <- c(a_accept, aux_a_new[[2]])
      } else {
        a_new <- aux_a_new
      }

    } else {
      a_new <- a_fix
      a_accept <- NULL
    }

    a <- a_new

    if(dev_verbose) {
      cat('...Done! \n')
    }

    # (g) Sampling "b"

    if(is.null(b_fix)){
      if(dev_verbose) {
        cat( 'Sampling "b":\n' )
      }
      aux_b_new <- sampling_b( n_sim_mh=1, b_ini=b,
                           a=a, d_0_b=d_0_b, d_1_b=d_1_b,
                           mu_star_n_r=mu_star_n_r,
                           max_it_time=max_it_time, burn_in=0, eta=eta,
                           accept_display=review_MH_chain, verbose=F )

      #plot(b_new)
      #hist(b_new)

      if (review_MH_chain) {
        b_new <- aux_b_new[[1]]
        b_accept <- c(b_accept, aux_b_new[[2]])
      } else {
        b_new <- aux_b_new
      }

      if(!all(b_new>-a)){stop('There is a problem sampling from "b", it should be >-a\nb=',b,"\n-a=",-a,sep="")}

    } else {
      b_new <- b_fix
      b_accept <- NULL
    }
    b <- b_new
    if(dev_verbose) {
      cat('...Done! \n')
    }

    # (h) Sampling "Z_ij"

    if(dev_verbose) {
      cat( 'Sampling "Z_ij":\n' )
    }
    Z_new <- get_latents( Y, data_class, mu_Z=mu_star[mu_star_map,], sigma_Z=sigma_Z, Z_old=Z)$Z
    #Z <- Z_new
    #browser()
    #Z <- get_latents( Y, data_class, mu_Z=mu_star[mu_star_map,], sigma_Z=sigma_Z, Z_old=Z )$Z
    if(dev_verbose) {
      cat('...Done! \n')
    }

    ### Storing simulation values for each iteration ###

    #mu_Z_sim[,,iter_i] <- mu_star[mu_star_map,]
    #sigma_Z_sim[,,iter_i] <- sigma_Z
    mu_star_map_sim[,paste("iter_",iter_i,sep="")] <- mu_star_map
    #mu_star_n_r_sim[[iter_i]] <- mu_star_n_r

    #browser()

    if(review_MH_chain) {
      Lambda_sim[iter_i,] <- diag(Lambda)
      Omega_sim[,,iter_i] <- Omega
      a_sim[iter_i] <- a
      b_sim[iter_i] <- b

    } else {
      Lambda_sim <- NULL
      Omega_sim[,,iter_i] <- NULL
      a_sim <- NULL
      b_sim <- NULL
    }

    if(!is.null(Z_sim_file)) {
      write.csv(Z,file=gsub(".csv",paste("_",formatC( iter_i, width = 2, flag = '0'),".csv",sep=""),Z_sim_file),quote=F,na="",row.names=F)
    }

    # save.image( file=paste(dir_out,'bayes_clustering.RData',sep='') )
    # load( file=paste(dir_out,'bayes_clustering.RData',sep='') )
    time_sim[iter_i+1] <- Sys.time()

    if(dev_verbose) {
      cat('\n   Finished iter_i=',iter_i,', \n     elapsed time:\n          ',format(diff(time_sim)[iter_i]),' this iteration\n          ',format(diff(time_sim[c(1,length(time_sim))])), ' total\n**----**',sep='')
    }

    if(!is.null(cluster_file)){
      if( is.element(iter_i,seq(0,n_iter,max(10,thin_k))) & iter_i>=(burn_in+thin_k) ) {
        iter_out <- seq(from=burn_in+1,to=iter_i,by=thin_k)
        bayes_cluster <- list(mu_star_map_sim=mu_star_map_sim[1+iter_out],

                              #mu_star_n_r_sim=mu_star_n_r_sim,
                              #mu_Z_sim=mu_Z_sim,
                              #sigma_Z_sim=sigma_Z_sim,
                              #mu_star_probs=mu_star_probs,

                              a_sim=a_sim[iter_out],
                              b_sim=b_sim[iter_out],
                              Lambda_sim=Lambda_sim[iter_out,],
                              Omega_sim=Omega_sim[,,iter_out],

                              Lambda_accept=Lambda_accept[iter_out,],
                              Omega_accept=Omega_accept[,,iter_out],
                              a_accept=a_accept[iter_out],
                              b_accept=b_accept[iter_out],

                              Lambda_accept_rate=apply(Lambda_accept[1:iter_i,],2,mean),
                              Omega_accept_rate=apply(Omega_accept[,,1:iter_i],c(1,2),mean),
                              a_accept_rate=mean(a_accept),
                              b_accept_rate=mean(b_accept)
        )
        save(bayes_cluster,file=cluster_file)
        rm(bayes_cluster); gc()
      }
    }

    if(diff(time_sim)[iter_i]>20*60) {
      cat("\n\n***** Clustering simulation ABORTED *****\n\n")
      if(dev_verbose) {
        cat('     Sys.time() = ',as.character(Sys.time()),sep='')
      }
      if( iter_i>=(burn_in+thin_k) ) {
        iter_out <- seq(from=burn_in+1,to=iter_i,by=thin_k)
        bayes_cluster <- list(mu_star_map_sim=mu_star_map_sim[1+iter_out],

                              #mu_star_n_r_sim=mu_star_n_r_sim,
                              #mu_Z_sim=mu_Z_sim,
                              #sigma_Z_sim=sigma_Z_sim,
                              #mu_star_probs=mu_star_probs,

                              a_sim=a_sim[iter_out],
                              b_sim=b_sim[iter_out],
                              Lambda_sim=Lambda_sim[iter_out,],
                              Omega_sim=Omega_sim[,,iter_out],

                              Lambda_accept=Lambda_accept[iter_out,],
                              Omega_accept=Omega_accept[,,iter_out],
                              a_accept=a_accept[iter_out],
                              b_accept=b_accept[iter_out],

                              Lambda_accept_rate=apply(Lambda_accept[1:iter_i,],2,mean),
                              Omega_accept_rate=apply(Omega_accept[,,1:iter_i],c(1,2),mean),
                              a_accept_rate=mean(a_accept),
                              b_accept_rate=mean(b_accept)
        )
        save(bayes_cluster,file=cluster_file)
        return(bayes_cluster)
      } else {
        cat('NO OUTPUT PRODUCED!!! \n the last iteration was still under the burn-in period!!! \n') ; return(NULL)
      }

    }

    rm(list=setdiff(ls(),lsa_iter))
    setTxtProgressBar(pb, iter_i)
    gc()
  }

  close(pb)

  if( iter_i>=(burn_in+thin_k) ) {
    iter_out <- seq(from=burn_in+1,to=iter_i,by=thin_k)
    bayes_cluster <- list(mu_star_map_sim=mu_star_map_sim[1+iter_out],

                          #mu_star_n_r_sim=mu_star_n_r_sim,
                          #mu_Z_sim=mu_Z_sim,
                          #sigma_Z_sim=sigma_Z_sim,
                          #mu_star_probs=mu_star_probs,

                          a_sim=a_sim[iter_out],
                          b_sim=b_sim[iter_out],
                          Lambda_sim=Lambda_sim[iter_out,],
                          Omega_sim=Omega_sim[,,iter_out],

                          Lambda_accept=Lambda_accept[iter_out,],
                          Omega_accept=Omega_accept[,,iter_out],
                          a_accept=a_accept[iter_out],
                          b_accept=b_accept[iter_out],

                          Lambda_accept_rate=apply(Lambda_accept[1:iter_i,],2,mean),
                          Omega_accept_rate=apply(Omega_accept[,,1:iter_i],c(1,2),mean),
                          a_accept_rate=mean(a_accept),
                          b_accept_rate=mean(b_accept)
    )

    #if(dev_verbose) {
    cat("\n*** Clustering estimation COMPLETED ***\n\n")
    #cat('     Sys.time() = ',as.character(Sys.time()),sep='')
    #cat('     Total time = ',as.character(Sys.time()-time_sim[1]),sep='')
    print(Sys.time()-time_sim[1])
    #}
    save(bayes_cluster,file=cluster_file)
    return(bayes_cluster)

  } else {
    cat('NO OUTPUT PRODUCED!!! \n the last iteration was still under the burn-in period!!! \n') ; return(NULL)
  }

}
