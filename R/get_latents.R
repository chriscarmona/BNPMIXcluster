
get_latents <- function( Y, data_class, mu_Z=NULL, sigma_Z=NULL, Z_old=NULL ) {

  ### This function simulates the latent variables for a specified Y ###
  #require(truncnorm)

  Y <- data.frame(as.matrix(Y))
  n <- nrow(Y)
  p <- ncol(Y)

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

  # number of variables by type
  n_c <- sum( is.element( data_class, var_classes[1] ) )
  n_o <- sum( is.element( data_class, var_classes[2] ) )
  n_m <- sum( is.element( data_class, var_classes[3] ) )

  p==n_c+n_o+n_m # TRUE

  # Sorting Y columns
  Y_new_order <- c( which(is.element(data_class,var_classes[1])),
                    which(is.element(data_class,var_classes[2])),
                    which(is.element(data_class,var_classes[3])) )
  data_class <- data_class[Y_new_order]
  Y <- Y[,Y_new_order]
  rm(Y_new_order)

  # changes the colnames of Y for simplicity and standarization
  colnames(Y) <- paste("var_",
                       c(rep(var_classes[1],n_c),rep(var_classes[2],n_o),rep(var_classes[3],n_m)),"_",
                       formatC( unlist(mapply(seq,1,c(n_c,n_o,n_m),length.out=c(n_c,n_o,n_m))) , width = 2, flag = '0'),
                       sep="")

  # Ordinal variables as factor
  if(n_o>0) {
    for (i in which( is.element( data_class, var_classes[2] ) ) ) {
      Y[,i] <- factor(Y[,i])
    }
  }

  # Categorical variables as factor #
  if(n_m>0) {
    n_m_l <- rep(NA,n_m) # vector with the number of latent variables that will be used for each categorical variable
    K_m <- rep(NA,n_m) # vector with the number of classes for each categorical variable
    for (i in 1:n_m ) {
      aux_i <- which( is.element( data_class, var_classes[3] ) )[i]
      Y[,aux_i] <- factor(Y[,aux_i])
      K_m[i] <- length(levels(Y[ ,aux_i])) # number of categories in that categorical variable
      n_m_l[i] <- K_m[i] - 1 # how many latents are needed for each categorical?
    }
  } else {n_m_l<-as.numeric(NULL)}

  # Compute the total number of latents that will be needed #
  #browser()
  n_q <- n_c+n_o+sum(n_m_l)

  ### Latent Variables ###
  Z <- matrix(data=NA,nrow=nrow(Y),ncol=n_q)

  # colnames of Z #
  colnames_Z <- NULL
  # colnames of Z: continuos #
  if( n_c > 0 ) {
    colnames_Z <- c( colnames_Z, paste("var_",
                                       rep(var_classes[1],n_c),"_",
                                       formatC( 1:n_c, width = 2, flag = '0'), sep="") )
  }
  # colnames of Z: ordinal #
  if( n_o > 0 ) {
    colnames_Z <- c( colnames_Z, paste("var_",
                                       rep(var_classes[2],n_o),"_",
                                       formatC( 1:n_o, width = 2, flag = '0'), sep="") )
  }
  # colnames of Z: categorical #
  if( n_m > 0 ) {
    for(i in 1:n_m) {
      colnames_Z <- c(colnames_Z, paste("var_",
                                        c( rep(var_classes[3],n_m_l[i] ) ),"_",
                                        formatC( i, width = 2, flag = '0'),"_",
                                        formatC( 1:n_m_l[i] , width = 2, flag = '0'),
                                        sep=""))
    }
  }
  colnames(Z) <- colnames_Z
  #browser()
  ##### Simulating Latent values #####

  if( is.null(mu_Z) ) {
    mu_Z <- matrix(0,nrow=n,ncol=n_q)
    colnames(mu_Z) <- colnames_Z
  }
  if( is.null(sigma_Z) ) {
    sigma_Z <- diag(1,nrow=n_q,ncol=n_q)
    colnames(sigma_Z) <- rownames(sigma_Z) <- colnames_Z
  }


  ###   Latent: Continuous   ###
  #browser()
  if(n_c>0) {
    # cat('Simulating latents for ',n_c,' Continuous variables   ')
    if( is.null(Z_old) ) {
      Z[,1:n_c] <- as.matrix( Y[,1:n_c] )
      # standardize the continuos variables
      #Z[,1:n_c] <- as.matrix( scale(Y[,1:n_c]) , nrow=nrow(Y), ncol=n_c )
    } else {
      Z[,1:n_c] <- Z_old[,1:n_c]
    }
  }


  ###   Latent: Ordinal   ###
  if(n_o>0) {
    #cat('Simulating latents for ',n_o,' Ordinal variables...')
    thres_o <- list() # list with the vector of thresholds for each ordinal variable
    K_o <- as.numeric(rep(NA,n_o)) # number of categories in each Ordinal variable

    # getting number of categories #

    for ( j in 1:n_o ) {
      # j<-1
      # cat(j,',')
      Y_ord_j <- which( is.element(data_class,var_classes[2])) # what columns in Y are ordinal?
      Y_ord_j <- Y_ord_j[j] # choose the ith ordinal

      # number of categories in that ordinal variable
      K_o[j] <- length(levels(Y[,Y_ord_j]))

      # obtain the thresholds dividing (-Inf,Inf) into K_o[i] intervals
      if(T) {
        # recommended thresholds of length 4 with a fixed variance of 1
        thres_o[[j]] <- c(-Inf,seq(from=-4*floor((K_o[j]-2)/2),to=4*ceiling((K_o[j]-2)/2),by=4),Inf)
      }
      if(F) {
        # equiprobable with respect to a standard normal distribution
        unif_int_j <- (0:K_o[j])/K_o[j] # unif_int_j are the equivalent thresholds for a uniform
        unif_int_j <- seq(0,1,by=1/K_o[j]) # unif_int_j are the equivalent thresholds for a uniform
        thres_o[[j]] <- qnorm(unif_int_j,0,1) # thres_o sets the actual thresholds for maping Y and Z
      }
      if(F) {
        # according to their distribution in the data
        unif_int_j <- c(0,cumsum(prop.table(table(Y[,Y_ord_j]))))
        thres_o[[j]] <- qnorm(unif_int_j,0,1) # thres_o sets the actual thresholds for maping Y and Z
      }
      # with the same lenght
      if(F) {
        thres_o[[j]] <- c(-Inf,0:(K_o[j]-2),Inf)
      }

    }

    # simulating ordinal latents #
    #if( !is.null(Z_old) ) {browser()}
    for(i in 1:n) {
      # i<-1
      for(j in 1:n_o) {
        # j<-1
        cat_ij <- match( Y[i,n_c+j],levels(Y[,n_c+j]) )

        if( is.null(Z_old) ) {

          if(F) {
            # simplified
            mu_Z_aux <- mu_Z[i,n_c+j]
            sigma_Z_aux <- sigma_Z[n_c+j,n_c+j]
            sigma_Z_aux <- sigma_Z_aux^(1/2)
          } else {
            # the point will come from a distribution such that
            # 95% of probability is within the interval
            # mu_Z_ij centered in its interval
            mu_Z_aux <- thres_o[[j]][cat_ij+c(0,1)]
            if(mu_Z_aux[1]==-Inf){mu_Z_aux[1]<-mu_Z_aux[2]-2*qnorm(0.95)}
            if(mu_Z_aux[2]==Inf){mu_Z_aux[2]<-mu_Z_aux[1]+2*qnorm(0.95)}
            mu_Z_aux <- mean(mu_Z_aux)
            mu_Z[i,n_c+j] <- mu_Z_aux

            # the variance is set such that it has 95% of being at that interval
            sigma_Z_aux <- thres_o[[j]][cat_ij+c(0,1)]
            if(any(abs(sigma_Z_aux)==Inf)) {
              sigma_Z_aux <- 1
            } else {
              sigma_Z_aux <- diff(sigma_Z_aux) / (2*qnorm(0.975))
            }

            #browser()
            # 95% of probability within the interval
            # pnorm(q=thres_o[[j]][cat_ij+c(0,1)],mean=mu_Z_aux,sd=sigma_Z_aux)
          }

        } else {
          #browser()
          mu_Z_aux <- mu_Z[i,n_c+j] + sigma_Z[(n_c+j),-(n_c+j)] %*% ginv(sigma_Z[-(n_c+j),-(n_c+j)]) %*% ( Z_old[i,-(n_c+j)] - mu_Z[i,-(n_c+j)] )
          sigma_Z_aux <- sigma_Z[n_c+j,n_c+j] + sigma_Z[(n_c+j),-(n_c+j)] %*% ginv(sigma_Z[-(n_c+j),-(n_c+j)]) %*% sigma_Z[-(n_c+j),(n_c+j)]
          sigma_Z_aux <- sigma_Z_aux^(1/2)
        }

        if(T) {
          Z[i,n_c+j] <- truncnorm::rtruncnorm(n=1,
                                   a=thres_o[[j]][cat_ij], b=thres_o[[j]][cat_ij+1],
                                   mean=mu_Z_aux,
                                   sd=sigma_Z_aux)

        }

        if( !(Z[i,n_c+j] > thres_o[[j]][cat_ij] & Z[i,n_c+j] < thres_o[[j]][cat_ij+1]) ){
          cat('\nError: There was a problem simulating an ordinal latent Z_ij, i=',i,', j=',n_c+j,'\n')
          stop('There was a problem simulating an ordinal latent Z_ij, i=',i,', j=',n_c+j,sep="")
        }
        #if( !(Z[i,n_c+j] > thres_o[[j]][cat_ij] & Z[i,n_c+j] < thres_o[[j]][cat_ij+1]) ){browser()}
      }

    }
    #cat('   done!\n')
  } else {
    thres_o<-as.numeric(NULL)
    K_o<-as.numeric(NULL)
  }

  ###   Latent: Categorical   ###

  if(n_m>0) {
    #cat('Simulating latents for ',n_m,' Categorical variables...')

    for(j in 1:n_m) {
      #cat(j,', ')
      # attention: j goes one for each categorical, not for each latent!

      # j<-2
      # defines dummies of Y_cat
      Y_j <- as.factor(Y[,n_c+n_o+j])
      Y_j <- relevel(Y_j,ref=length(levels(Y_j))) # assign the reference level as the last level
      Y_j <- model.matrix( ~.,data = data.frame( Y_j )) # tranform the categorical to dummy variables
      Y_j <- as.matrix( Y_j[,-1] )

      #head(Y_j)
      #head(Y[,n_c+n_o+j])

      for ( i in 1:n ) {
        # i<-2
        # simulating one by one
        if( any( Y_j[i,]>0 ) ) { # observation i is different from the last class, at least one latent >0 and the max corresponds to the value of the categorical

          # first, simulates the maximum #
          j_z_max <- n_c + n_o + diffinv(n_m_l)[j] + which(Y_j[i,]==1) # column of z corresponding with the maximum

          if( is.null(Z_old) ) {
            mu_Z_aux <- mu_Z[i,j_z_max]
            sigma_Z_aux <- sigma_Z[j_z_max,j_z_max]
            sigma_Z_aux <- sigma_Z_aux^(1/2)
          } else {
            mu_Z_aux <- mu_Z[i,j_z_max] + sigma_Z[j_z_max,-j_z_max] %*% ginv(sigma_Z[-j_z_max,-j_z_max]) %*% ( Z_old[i,-j_z_max] - mu_Z[i,-j_z_max] )
            sigma_Z_aux <- sigma_Z[j_z_max,j_z_max] + sigma_Z[j_z_max,-j_z_max] %*% ginv(sigma_Z[-j_z_max,-j_z_max]) %*% sigma_Z[-j_z_max,j_z_max]
            sigma_Z_aux <- sigma_Z_aux^(1/2)
          }

          sim_int <- c(pnorm(0,mean=mu_Z_aux,sd=sigma_Z_aux),1) # simulates a positive gaussian number for the maximum
          Z[i, j_z_max ] <- qnorm( p=runif(n=1, min=sim_int[1], max=sim_int[2]) , mean=mu_Z_aux, sd=sigma_Z_aux)

          # then, simulates the rest #
          for( s in (1:ncol(Y_j))[-which(Y_j[i,]==1)] ) {
            j_z_s <- n_c + n_o + diffinv(n_m_l)[j] + s # column of z corresponding with this value
            if( is.null(Z_old) ) {
              mu_Z_aux <- mu_Z[i,j_z_s]
              sigma_Z_aux <- sigma_Z[j_z_s,j_z_s]
              sigma_Z_aux <- sigma_Z_aux^(1/2)
            } else {
              mu_Z_aux <- mu_Z[i,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% ginv(sigma_Z[-j_z_s,-j_z_s]) %*% ( Z_old[i,-j_z_s] - mu_Z[i,-j_z_s] )
              sigma_Z_aux <- sigma_Z[j_z_s,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% ginv(sigma_Z[-j_z_s,-j_z_s]) %*% sigma_Z[-j_z_s,j_z_s]
              sigma_Z_aux <- sigma_Z_aux^(1/2)
            }

            sim_int <- c(0,pnorm( Z[i, j_z_max ] ,mean=mu_Z_aux,sd=sigma_Z_aux)) # simulates a gaussian number lower than the maximum
            Z[i, j_z_s ] <- qnorm( p=runif(n=1, min=sim_int[1], max=sim_int[2]) , mean=mu_Z_aux, sd=sigma_Z_aux)
          }

          # checking consistency #
          if( which(Y_j[i,]==1) != which(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)]==max(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)])) ) {
            cat('\nError: There is a problem generating a categorical variable i=',i,' j=',j,'\n')
            stop('There is a problem generating a categorical variable i=',i,' j=',j,sep='')
          }

        } else { # observation i is equal to the last class, all values of latents have to be negative
          for( s in 1:ncol(Y_j) ) {
            j_z_s <- n_c + n_o + diffinv(n_m_l)[j] + s # column of z corresponding with this value

            if( is.null(Z_old) ) {
              mu_Z_aux <- mu_Z[i,j_z_s]
              sigma_Z_aux <- sigma_Z[j_z_s,j_z_s]
              sigma_Z_aux <- sigma_Z_aux^(1/2)
            } else {
              mu_Z_aux <- mu_Z[i,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% ginv(sigma_Z[-j_z_s,-j_z_s]) %*% ( Z_old[i,-j_z_s] - mu_Z[i,-j_z_s] )
              sigma_Z_aux <- sigma_Z[j_z_s,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% ginv(sigma_Z[-j_z_s,-j_z_s]) %*% sigma_Z[-j_z_s,j_z_s]
              sigma_Z_aux <- sigma_Z_aux^(1/2)
            }

            sim_int <- c(0,pnorm(0,mean=mu_Z_aux,sd=sigma_Z_aux)) # simulates negative gaussian numbers
            Z[i, j_z_s ] <- qnorm( p=runif(n=1, min=sim_int[1], max=sim_int[2]) , mean=mu_Z_aux, sd=sigma_Z_aux)
          }
          # checking consistency #
          all(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)]<0)
        }
      }
    }
    #cat('...done!\n')
  }

  dim(Z)[1]==n # TRUE
  dim(Z)[2]==n_q # TRUE

  return( list(Z=Z,
               n_c=n_c,n_o=n_o,n_m_l=n_m_l,
               n_q=n_q,
               thres_o=thres_o,K_o=K_o) )
}
