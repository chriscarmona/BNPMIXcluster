#'
#' @title
#'     Simulation of latent variables \eqn{Z} in the \emph{mixdpclust} model
#'
#' @description
#'     Simulates values for latent variables \eqn{Z=(Z_1,...,Z_q)} according to the specification in the \emph{mixdpclust} model.
#'
#' @param Y Matrix or data frame containing the observed data.
#'
#' @param var_type Character vector that indicates the type of variable in each column of Y. Three possible types:
#' \itemize{
#'   \item "\strong{c}" for continuous variables.
#'   \item "\strong{o}" for ordinal variables (ordered categorical).
#'   \item "\strong{m}" for nominal variables (non-ordered categorical).
#' }
#'
#' @param mu_Z.ini an optional vector with the expected values \eqn{\mu_Z} of the latent variables.
#' @param sigma_Z.ini an optional matrix with the covariance matrix \eqn{\Sigma_Z} of the latent variables.
#' @param Z.ini an optional matrix with initial values for the latent variables that will be simulated.
#'
#' @importFrom stats diffinv model.matrix pnorm qnorm relevel runif
#' @importFrom MASS ginv
#'
#' @details
#'     For each variable in the \code{Y} data frame, an associated continuos latent variable is generated.
#'     if \code{var_type} is continuos, the corresponding \eqn{Z} will keep the original values of \eqn{Y}.
#'     If \eqn{Y} is categorical, the function will scan the unique values of \eqn{Y} and generate continuous latent variables accordingly.
#'
#' @keywords internal

get_latents <- function( Y,
                         var_type,
                         mu_Z.ini=NULL,
                         sigma_Z.ini=NULL,
                         Z.ini=NULL ) {

  ### This function simulates the latent variables for a specified Y ###
  #require(truncnorm)

  Y <- data.frame(as.matrix(Y))
  n.obs <- nrow(Y)
  p <- ncol(Y)

  # possible variable classes that are allowed
  var_type_all <- c("c","o","m")

  # checking input consistency
  if( length(var_type) != ncol(Y) ) {
    cat('\nError: The number of columns in "Y" have to be equal to the lenght of vector "var_type"\n')
    stop('The number of columns in "Y" have to be equal to the lenght of vector "var_type"')
  }
  if( any(!is.element(var_type,var_type_all)) ) {
    cat('\nError: Elements in "var_type" have to be one of ',paste(var_type_all,collapse = ","),'\n')
    stop('Elements in "var_type" have to be one of ',paste(var_type_all,collapse = ","))
  }

  # number of variables by type
  n_c <- sum( is.element( var_type, var_type_all[1] ) )
  n_o <- sum( is.element( var_type, var_type_all[2] ) )
  n_m <- sum( is.element( var_type, var_type_all[3] ) )

  p==n_c+n_o+n_m # TRUE

  # Sorting Y columns
  Y_new_order <- c( which(is.element(var_type,var_type_all[1])),
                    which(is.element(var_type,var_type_all[2])),
                    which(is.element(var_type,var_type_all[3])) )
  var_type <- var_type[Y_new_order]
  Y <- Y[,Y_new_order]
  rm(Y_new_order)

  # changes the colnames of Y for simplicity and standarization
  colnames(Y) <- paste("var_",
                       c(rep(var_type_all[1],n_c),rep(var_type_all[2],n_o),rep(var_type_all[3],n_m)),"_",
                       formatC( unlist(mapply(seq,1,c(n_c,n_o,n_m),length.out=c(n_c,n_o,n_m))) , width = 2, flag = '0'),
                       sep="")

  # Transform ordinal variables as factor
  if(n_o>0) {
    for (i in which( is.element( var_type, var_type_all[2] ) ) ) {
      Y[,i] <- factor(Y[,i])
    }
  }

  if(n_m>0) {

    # vector with the number of latent variables that will be used for each categorical variable
    n_m_l <- rep(NA,n_m)

    # vector with the number of classes for each categorical variable
    K_m <- rep(NA,n_m)
    for (i in 1:n_m ) {
      # which columns are nominal categorical?
      aux_i <- which( is.element( var_type, var_type_all[3] ) )[i]

      # Transfrom categorical variables as factor #
      Y[,aux_i] <- factor(Y[,aux_i])

      # number of categories in that categorical variable
      K_m[i] <- length(levels(Y[ ,aux_i]))

      # how many latents are needed for each categorical?
      n_m_l[i] <- K_m[i] - 1
    }
    rm(i,aux_i)

  } else {
    n_m_l<-as.numeric(NULL)
  }

  # Compute the total number of latents that will be needed #
  n_q <- n_c+n_o+sum(n_m_l)

  ### Latent Variables ###
  Z <- matrix(data=NA,nrow=n.obs,ncol=n_q)

  # colnames of Z #
  colnames_Z <- NULL

  # colnames of Z: continuos #
  if( n_c > 0 ) {
    colnames_Z <- c( colnames_Z, paste("var_",
                                       rep(var_type_all[1],n_c),"_",
                                       formatC( 1:n_c, width = 2, flag = '0'), sep="") )
  }
  # colnames of Z: ordinal #
  if( n_o > 0 ) {
    colnames_Z <- c( colnames_Z, paste("var_",
                                       rep(var_type_all[2],n_o),"_",
                                       formatC( 1:n_o, width = 2, flag = '0'), sep="") )
  }
  # colnames of Z: categorical #
  if( n_m > 0 ) {
    for(i in 1:n_m) {
      colnames_Z <- c(colnames_Z, paste("var_",
                                        c( rep(var_type_all[3],n_m_l[i] ) ),"_",
                                        formatC( i, width = 2, flag = '0'),"_",
                                        formatC( 1:n_m_l[i] , width = 2, flag = '0'),
                                        sep=""))
    }
  }
  colnames(Z) <- colnames_Z

  ### Mean and Variance of Z ###

  if( is.null(mu_Z.ini) ) {
    mu_Z <- matrix(0,nrow=n.obs,ncol=n_q)
  } else {
    mu_Z <- mu_Z.ini
  }
  colnames(mu_Z) <- colnames_Z

  if( is.null(sigma_Z.ini) ) {
    sigma_Z <- diag(1,nrow=n_q,ncol=n_q)
  } else {
    sigma_Z <- sigma_Z.ini
  }
  colnames(sigma_Z) <- rownames(sigma_Z) <- colnames_Z

  ##### Simulating Latent values #####

  ###   Latent: Continuous   ###
  # Continuous variables will remain the same values #
  if(n_c>0) {
    # cat('Simulating latents for ',n_c,' Continuous variables   ')
    if( is.null(Z.ini) ) {
      # Z[,1:n_c] <- as.matrix( Y[,1:n_c] )
      # standardize the continuos variables
      Z[,1:n_c] <- as.matrix( scale(Y[,1:n_c]) , nrow=n.obs, ncol=n_c )
    } else {
      Z[,1:n_c] <- Z.ini[,1:n_c]
    }
  }


  ###   Latent: Ordinal   ###
  if(n_o>0) {

    # list with the vector of thresholds for each ordinal variable
    thres_o <- list()

    # number of categories in each Ordinal variable
    K_o <- as.numeric(rep(NA,n_o))

    # getting number of categories #
    for ( j in 1:n_o ) {
      # what columns in Y are ordinal?
      Y_ord_j <- which( is.element(var_type,var_type_all[2]))
      # choose the ith ordinal
      Y_ord_j <- Y_ord_j[j]

      # number of categories in that ordinal variable
      K_o[j] <- length(levels(Y[,Y_ord_j]))

      # obtain the thresholds dividing (-Inf,Inf) into K_o[i] intervals

      # recommended thresholds of length 4 with a fixed variance of 1
      thres_o[[j]] <- c(-Inf,seq(from=-4*floor((K_o[j]-2)/2),to=4*ceiling((K_o[j]-2)/2),by=4),Inf)
    }

    # simulating ordinal latents #
    for(j in 1:n_o) {
      # j goes one for each ordinal in Y,  which is the same as each ordinal latent in Z!

      # j<-1
      for(i in 1:n.obs) {
        # i<-1

        # corresponding category of the choosen Y[i,n_c+j]
        cat_ij <- match( Y[i,n_c+j],levels(Y[,n_c+j]) )

        if( is.null(Z.ini) ) {
          # if there is no given initial value for Z[i,n_c+j]
          if(T) {
            # the simulated ordinal latent Z[i,n_c+j]
            # come from a distribution with mean mu_Z[i,n_c+j] and variance sigma_Z[n_c+j,n_c+j]
            mu_Z_ij <- mu_Z[i,n_c+j]
            sigma_Z_ij <- sigma_Z[n_c+j,n_c+j]
            sigma_Z_ij <- sigma_Z_ij^(1/2)
          } else {
            # the simulated ordinal latent Z[i,n_c+j]
            # come from a distribution such that
            # 95% of probability is within the given interval
            # mu_Z_ij centered in its interval
            mu_Z_ij <- thres_o[[j]][cat_ij+c(0,1)]
            if(mu_Z_ij[1]==-Inf){mu_Z_ij[1]<-mu_Z_ij[2]-2*qnorm(0.95)}
            if(mu_Z_ij[2]==Inf){mu_Z_ij[2]<-mu_Z_ij[1]+2*qnorm(0.95)}
            mu_Z_ij <- mean(mu_Z_ij)
            mu_Z[i,n_c+j] <- mu_Z_ij

            # the variance is set such that it has 95% of being at that interval
            sigma_Z_ij <- thres_o[[j]][cat_ij+c(0,1)]
            if(any(abs(sigma_Z_ij)==Inf)) {
              sigma_Z_ij <- 1
            } else {
              sigma_Z_ij <- diff(sigma_Z_ij) / (2*qnorm(0.975))
            }

            # 95% of probability within the interval
            # pnorm(q=thres_o[[j]][cat_ij+c(0,1)],mean=mu_Z_ij,sd=sigma_Z_ij)
          }

        } else {
          # if an initial value for Z was given
          mu_Z_ij <- mu_Z[i,n_c+j] + sigma_Z[(n_c+j),-(n_c+j)] %*% MASS::ginv(sigma_Z[-(n_c+j),-(n_c+j)]) %*% ( Z.ini[i,-(n_c+j)] - mu_Z[i,-(n_c+j)] )
          sigma_Z_ij <- sigma_Z[n_c+j,n_c+j] + sigma_Z[(n_c+j),-(n_c+j)] %*% MASS::ginv(sigma_Z[-(n_c+j),-(n_c+j)]) %*% sigma_Z[-(n_c+j),(n_c+j)]
          sigma_Z_ij <- sigma_Z_ij^(1/2)
        }

        Z[i,n_c+j] <- truncnorm::rtruncnorm(n=1,
                                            a=thres_o[[j]][cat_ij], b=thres_o[[j]][cat_ij+1],
                                            mean=mu_Z_ij,
                                            sd=sigma_Z_ij)

        if( !(Z[i,n_c+j] > thres_o[[j]][cat_ij] & Z[i,n_c+j] < thres_o[[j]][cat_ij+1]) ){
          cat('\nError: There was a problem simulating an ordinal latent Z_ij, i=',i,', j=',n_c+j,'\n')
          stop('There was a problem simulating an ordinal latent Z_ij, i=',i,', j=',n_c+j,sep="")
        }
      }
    }
    rm(i,j,mu_Z_ij,sigma_Z_ij)

  } else {
    thres_o<-as.numeric(NULL)
    K_o<-as.numeric(NULL)
  }

  ###   Latent: Categorical   ###

  if(n_m>0) {

    for(j in 1:n_m) {
      # j goes one for each nominal variable in Y, not for each nominal latent in Z!

      # defines dummies of Y_cat
      Y_j <- as.factor(Y[,n_c+n_o+j])

      # assign the reference level as the last level
      Y_j <- stats::relevel(Y_j,ref=length(levels(Y_j)))

      # tranform the categorical to dummy variables
      Y_j <- stats::model.matrix( ~.,data = data.frame( Y_j ))
      Y_j <- as.matrix( Y_j[,-1] )

      #head(Y_j)
      #head(Y[,n_c+n_o+j])

      for ( i in 1:n.obs ) {
        # i<-2
        # simulating one by one
        if( any( Y_j[i,]>0 ) ) {
          # observation Y[i,j] IS NOT in the last class

          # in this case:
          # at least one latent >0
          # and the maximum value of the associated Z's corresponds to the value of the class in Y

          ### first, simulates the maximum ###

          # column of Z[i,] corresponding with the category of Y[i,]
          # it will take the maximum value of the corresponding latent
          j_z_max <- n_c + n_o + diffinv(n_m_l)[j] + which(Y_j[i,]==1)

          if( is.null(Z.ini) ) {
            mu_Z_ij <- mu_Z[i,j_z_max]
            sigma_Z_ij <- sigma_Z[j_z_max,j_z_max]
            sigma_Z_ij <- sigma_Z_ij^(1/2)
          } else {
            mu_Z_ij <- mu_Z[i,j_z_max] + sigma_Z[j_z_max,-j_z_max] %*% MASS::ginv(sigma_Z[-j_z_max,-j_z_max]) %*% ( Z.ini[i,-j_z_max] - mu_Z[i,-j_z_max] )
            sigma_Z_ij <- sigma_Z[j_z_max,j_z_max] + sigma_Z[j_z_max,-j_z_max] %*% MASS::ginv(sigma_Z[-j_z_max,-j_z_max]) %*% sigma_Z[-j_z_max,j_z_max]
            sigma_Z_ij <- sigma_Z_ij^(1/2)
          }

          # simulates a positive gaussian number for the maximum
          Z[i, j_z_max ] <- truncnorm::rtruncnorm(n=1,
                                                a=0, b=Inf,
                                                mean=mu_Z_ij,
                                                sd=sigma_Z_ij)

          # then, simulates the rest #
          for( s in (1:ncol(Y_j))[-which(Y_j[i,]==1)] ) {
            j_z_s <- n_c + n_o + diffinv(n_m_l)[j] + s # column of z corresponding with this value
            if( is.null(Z.ini) ) {
              mu_Z_ij <- mu_Z[i,j_z_s]
              sigma_Z_ij <- sigma_Z[j_z_s,j_z_s]
              sigma_Z_ij <- sigma_Z_ij^(1/2)
            } else {
              mu_Z_ij <- mu_Z[i,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% MASS::ginv(sigma_Z[-j_z_s,-j_z_s]) %*% ( Z.ini[i,-j_z_s] - mu_Z[i,-j_z_s] )
              sigma_Z_ij <- sigma_Z[j_z_s,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% MASS::ginv(sigma_Z[-j_z_s,-j_z_s]) %*% sigma_Z[-j_z_s,j_z_s]
              sigma_Z_ij <- sigma_Z_ij^(1/2)
            }

            # simulates a gaussian number lower than the maximum
            Z[i, j_z_s ] <- truncnorm::rtruncnorm(n=1,
                                                  a=-Inf, b=Z[i, j_z_max ],
                                                  mean=mu_Z_ij,
                                                  sd=sigma_Z_ij)

          }

          # checking consistency #
          if( which(Y_j[i,]==1) != which(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)]==max(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)])) ) {
            cat('\nError: There is a problem generating a categorical variable i=',i,' j=',j,'\n')
            stop('There is a problem generating a categorical variable i=',i,' j=',j,sep='')
          }


        } else {
          # observation Y[i,j] IS in the last class

          # in this case:
          # all values of latents have to be negative

          for( s in 1:ncol(Y_j) ) {
            j_z_s <- n_c + n_o + diffinv(n_m_l)[j] + s # column of z corresponding with this value

            if( is.null(Z.ini) ) {
              mu_Z_ij <- mu_Z[i,j_z_s]
              sigma_Z_ij <- sigma_Z[j_z_s,j_z_s]
              sigma_Z_ij <- sigma_Z_ij^(1/2)
            } else {
              mu_Z_ij <- mu_Z[i,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% MASS::ginv(sigma_Z[-j_z_s,-j_z_s]) %*% ( Z.ini[i,-j_z_s] - mu_Z[i,-j_z_s] )
              sigma_Z_ij <- sigma_Z[j_z_s,j_z_s] + sigma_Z[j_z_s,-j_z_s] %*% MASS::ginv(sigma_Z[-j_z_s,-j_z_s]) %*% sigma_Z[-j_z_s,j_z_s]
              sigma_Z_ij <- sigma_Z_ij^(1/2)
            }

            # simulates negative gaussian numbers #
            Z[i, j_z_s ] <- truncnorm::rtruncnorm(n=1,
                                                  a=-Inf, b=0,
                                                  mean=mu_Z_ij,
                                                  sd=sigma_Z_ij)
          }

          # checking consistency of simulated values #
          if(!all(Z[i,n_c + n_o + diffinv(n_m_l)[j]+1:ncol(Y_j)]<0)) {
            stop("There is a problem with the simulation of nominal latents in Z")
          }

        }
      }
    }

  }

  if(dim(Z)[1]!=n.obs) {
    stop("There is a problem with the simulation of latents Z")
  }
  if(dim(Z)[2]!=n_q) {
    stop("There is a problem with the simulation of latents Z")
  }

  latents_data <- list( Z=Z,
                  n_c=n_c,
                  n_o=n_o,
                  n_m_l=n_m_l,
                  n_q=n_q,
                  thres_o=thres_o,
                  K_o=K_o )

  return( latents_data )
}
