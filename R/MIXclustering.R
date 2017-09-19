#' @title
#'     Bayesian Nonparametric Model for Clustering with Mixed Scale Variables
#'
#' @description
#'     \code{MIXclustering} is used to perform cluster analisis of individuals using a Bayesian nonparametric mixture model that jointly models mixed scale data and accommodates for different sampling probabilities. The model is described in Carmona, C., Nieto-Barajas, L. E., Canale, A. (2016).
#'
#' @param Y Matrix or data frame containing the data to be clustered.
#'
#' @param var_type Character vector that indicates the type of variable in each column of x. Three possible types:
#' \itemize{
#'   \item "\strong{c}" for continuous variables. It is assumed to be Gaussian-shaped.
#'   \item "\strong{o}" for ordinal variables (binary and ordered categorical).
#'   \item "\strong{m}" for nominal variables (non-ordered categorical).
#' }
#'
#' @param sampling_prob vector with the sampling probabilities \eqn{\pi_i} for each individual in case that the data come from a complex survey sample. By default \eqn{\pi_i=1}.
#' @param expansion_f vector with the expansion factors, the reciprocal of the sampling probabilities, \eqn{w_i = 1/\pi_i}. If both \code{sampling_prob} and \code{expansion_f} are specified, preference is given to \code{sampling_prob}.
#'
#' @param n_iter_out Number of effective iterations in the MCMC procedure for clustering.
#' @param n_burn Number of iterations discarded as part of the burn-in period at the beginning MCMC procedure.
#' @param n_thin Number of iterations discarded for thining the chain (reducing the autocorrelation). We keep 1 of every n_thin iterations.
#'
#' @param a_fix A numeric value to set the parameter \eqn{a} in the model. If \code{NULL} (default), the parameter \eqn{a} is assigned a prior distribution. See \code{details}.
#' @param alpha Hyperparameter in the prior distribution of \eqn{a}. See \code{details}.
#' @param d_0_a Hyperparameter in the prior distribution of \eqn{a}. See \code{details}.
#' @param d_1_a Hyperparameter in the prior distribution of \eqn{a}. See \code{details}.
#'
#' @param b_fix A numeric value to set the parameter \eqn{b} in the model. If \code{NULL} (default), the parameter \eqn{b} is assigned a prior distribution. See \code{details}.
#' @param d_0_b Hyperparameter in the prior distribution of \eqn{b}. See \code{details}.
#' @param d_1_b Hyperparameter in the prior distribution of \eqn{b}. See \code{details}.
#' @param eta Tuning parameter controlling the proposal in the \emph{Metropolis-Hastings} step for \eqn{b}.
#'
#' @param d_0_z Hyperparameter in the prior distribution of the variance for the latent variables. See \code{details}.
#' @param d_1_z Hyperparameter in the prior distribution of the variance for the latent variables. See \code{details}.
#' @param kappa Tuning parameter controlling the proposal in the \emph{Metropolis-Hastings} step for the variance of latent variables.
#'
#' @param delta Tuning parameter controlling the proposal in the \emph{Metropolis-Hastings} step for the correlation of latent variables.
#'
#' @param d_0_mu Hyperparameter in the prior distribution of the variance of the location in each cluster. See \code{details}.
#' @param d_1_mu Hyperparameter in the prior distribution of the variance of the location in each cluster. See \code{details}.
#'
#' @param log_file Specifies a file to save the details with the execution time and the parameters used.
#'
#' @details
#'
#' The model consists on a bayesian non-parametric approach for clustering that is capable to combine different types of variables through the usage of associated continuous latent variables. The clustering mechanism is based on a location mixture model with a Poisson-Dirichlet (\eqn{PD}) process prior on the location parameters \eqn{\mu_i ; i=1,...,n} of the associated latent variables.
#'
#' Computational inference about the cluster allocation and the posterior distribution of the parameters are performed using MCMC simulations.
#'
#' A full description of the model is in the article Carmona et al. (2016) (preprint: \url{http://arxiv.org/abs/1612.00083}). See \code{Reference}.
#'
#' The model consider an individual \eqn{y_i} that is characterized by a multivariate response of dimension \eqn{p}, i.e., \eqn{y_i=(y_{i,1},...,y_{i,p})}. The total number of variables \eqn{p} is divided into \eqn{c} continuous variables, \eqn{o} ordinal variables, and \eqn{m} nominal variables such that \eqn{p=c+o+m}.
#'
#' For the continuous variables, it is convenient that the variables have a real support. The user may have transformed the original values before using the function \code{MIXclustering}.
#'
#' For each response \eqn{y_i=(y_{i,1},...,y_{i,p})} (of dimension \eqn{p}) a corresponding latent vector \eqn{z_i=(z_{i,1},...,z_{i,q})} (of dimension \eqn{q}) is created, according to the following:
#' \itemize{
#'    \item For each continuous variable \eqn{y_{i,j} ; j=1,...,c} the algorithm uses a latent with the same values \eqn{z_{i,j}=y_{i,j}}.
#'    \item For each ordinal variable \eqn{y_{i,j} , j = c+1,...,c+o}, with \eqn{K_j} different ordered values, the algorithm creates one latent \eqn{z_{i,j}}, that allows to map the categories into continuous values divided by thresholds. For example, for a binary \eqn{y_j}, we have \eqn{y_j=0} iff \eqn{z_j<0} and \eqn{y_j=1} iff \eqn{z_j>0}
#'    \item For each nominal variable \eqn{y_{i,j} , j = c+o+1,...,c+o+m}, with \eqn{L_j} categories, the algorithm require \eqn{L_j-1} latent variables, whose relative order is consistent with the observed category.
#'    }
#'
#' The data may come from a complex survey sample where each individual \eqn{y_i} has known sampling probability \eqn{\pi_i, i = 1,...,n}. The reciprocal of these sampling probabilities, \eqn{w_i=1/\pi_i}, are called expansion factors or sampling design weights.
#'
#' The joint model for the latent vector is therefore:
#'
#' \eqn{(z_i | \mu_i,\Sigma)} ~ \eqn{N_q(\mu_i, \pi_i \Sigma )}
#'
#' (Note: the final model in Carmona et al. (2016) has variance \eqn{\kappa \pi_i \Sigma}. This value of \eqn{\kappa} can be used in the package through a transformed sampling probability vector \eqn{\pi_i^*=\kappa\pi_i})
#'
#' The clustering model will be based in an appropriate choice of the prior distribution on the \eqn{\mu_i}'s. A clustering of the \eqn{\mu_i}'s will induce a clustering of the \eqn{y_i}'s. Our prior on the \eqn{\mu_i}'s will be:
#'
#' \eqn{\mu_i | G}~\eqn{G}, iid for \eqn{i=1,...,n}
#'
#' Where \eqn{G}~\eqn{PD(a,b,G_0)} is a Poisson-Dirichlet process with parameters \eqn{a \in [0,1)}, \eqn{b>-a} and centring measure \eqn{G_0}. The Dirichlet and the normalized stable processes arise when \eqn{a=0} and when \eqn{b=0}, respectively.
#'
#' In consequence, this choice of prior implies that the \eqn{\mu_i}'s are exchangeable with marginal distribution \eqn{\mu_i}~\eqn{G_0} for all \eqn{i=1,...,n}.
#'
#' In our case, \eqn{G(\mu)=N(0,\Sigma_\mu)}, where \eqn{\Sigma_\mu = diag( \sigma^2_{\mu 1} ,...,\sigma^2_{\mu q} )}.
#'
#' The parameters \eqn{a} and \eqn{b} in the model define the \eqn{PD} process and therefore control the number of groups. These parameters can be fixed, resulting in a larger/smaller number of groups if assigned a larger/smaller value, respectively.
#'
#' There are 9 hyperparameters in the function that also characterize the prior distributions in the model:
#' \itemize{
#'     \item f(a) = alpha * I(a=0) + (1-alpha) * dbeta( a | d_0_a , d_0_a )
#'     \item f(b | a) = dgamma( b + a | d_0_b , d_1_b )
#'     \item sigma^2 ~ inverse-gamma( d_0_z , d_1_z)
#'     \item sigma^2_mu ~ inverse-gamma( d_0_mu , d_1_mu )
#' }
#'
#' The definition of these values also affect the number of resulting clusters since they affect the variance implied in the model.
#'
#' For example, increasing the values of \code{d_1_a} and \code{d_1_b} reduce the number of groups.
#'
#' Finally, the function parameters \code{eta}, \code{kappa}, \code{delta} are tuning parameters that control the acceptance rate in the random-walk MH steps of the new proposed values for the parameters \eqn{b}, \eqn{\Lambda_{j,j}} (variance of latents) and \eqn{\Omega_{i,j}} (correlation of latents). These parameters are not recomended to be changed (used in the internal functions: \code{sampling_b} , \code{sampling_Lambda_jj} , \code{sampling_Omega_ij}).
#'
#' @return
#' \code{MIXclustering} returns a S3 object of class "\code{MIXcluster}".
#'
#' The generic methods \code{\link{summary}} and \code{\link{plot}} are defined for this class.
#'
#' An object of class "MIXcluster" is a list containing the following components:
#' \describe{
#'     \item{\code{cluster}}{vector with the cluster allocation for each row in the data. It corresponds to the iteration which is Closest-To-Average (CTA) arrangement.}
#'     \item{\code{cluster_heterogeneity}}{Heterogeneity Measure (HM) for the cluster in the previous point. The HM measure is discussed in section 4 of Carmona et al. (2017).}
#'     \item{\code{Y.cluster.summary}}{a summary of the data divided by the allocation in \code{$cluster}.}
#'     \item{\code{Y.var_type}}{vector with the variable types in the data.}
#'     \item{\code{Y.na}}{vector specifying the rows with missing values.}
#'     \item{\code{Y.n}}{number of rows in the data.}
#'     \item{\code{Y.p}}{number of variables in the data.}
#'     \item{\code{MC.clusters}}{matrix with the cluster allocation for each row in the data. Each column corresponds to an effective iteration in the MCMC simulation of the model (after discarding burn-in and thinning iterations).}
#'     \item{\code{MC.clusters_heterogeneity}}{Heterogeneity Measure (HM) for all the clusters returned in \code{MC.clusters}.}
#'     \item{\code{cluster.matrix.avg}}{average similarity matrix of size \eqn{n} by \eqn{n}.}
#'     \item{\code{MC.values}}{a list with the simulated values of the chains for the parameters \eqn{a},\eqn{b},\eqn{\Lambda},\eqn{\Omega}.}
#'     \item{\code{MC.accept.rate}}{a named vector with the acceptance rates for each parameter. It includes iterations that are discarded in the burn-in period and thinning.}
#'     \item{\code{call}}{the matched call.}
#' }
#'
#' @seealso
#' \code{\link{summary.MIXcluster}} for a summary of the clustering results, \code{\link{plot.MIXcluster}} for graphical representation of results.
#'
#' @examples
#' 
#' ##############################
#' #     Simulation study 1     #
#' #    Carmona et al. (2017)   #
#' ##############################
#' 
#' # Data and parameters are discussed in section 5.1 of Carmona et al. (2017) #
#' 
#' # Set seed for reproducibility #
#' set.seed(0) 
#' 
#' 
#' # Specification of data Y #
#' help(Y_ex_5_1)
#' # Choose 1, 2, or 3 #
#' ex_i <- 1
#' 
#' # specification of parameters #
#' help(meta_param_ex)
#' # Choose "a", "b" or "c" #
#' param_j <- "c"
#' 
#' var_type_Y_ex_5_1 <- list( c("c","c","c"),
#'                            c("o","o"),
#'                            c("o","o","o","c") )
#' \dontrun{
#' cluster_ex <- MIXclustering( Y = as.matrix(Y_ex_5_1[[ ex_i ]]),
#'                              var_type=var_type_Y_ex_5_1[[ ex_i ]],
#'                              
#'                              n_iter_out=1500,
#'                              n_burn=200,
#'                              n_thin=3,
#'                              
#'                              alpha = meta_param_ex[ param_j, "alpha" ],
#'                              d_0_a = meta_param_ex[ param_j, "d_0_a"],
#'                              d_1_a = meta_param_ex[ param_j, "d_1_a" ],
#'                              d_0_b = meta_param_ex[ param_j, "d_0_b" ],
#'                              d_1_b = meta_param_ex[ param_j, "d_1_b" ],
#'                              eta = meta_param_ex[ param_j, "eta" ],
#'                              kappa = meta_param_ex[ param_j, "kappa" ],
#'                              delta = meta_param_ex[ param_j, "delta" ],
#'                              
#'                              d_0_z = meta_param_ex[ param_j, "d_0_z" ],
#'                              d_1_z = meta_param_ex[ param_j, "d_1_z" ],
#'                              d_0_mu = meta_param_ex[ param_j, "d_0_mu" ],
#'                              d_1_mu = meta_param_ex[ param_j, "d_1_mu" ] )
#' # Summary of clustering results
#' summary(cluster_ex)
#'
#' # Representation of clustering results
#' plot(cluster_ex,type="heatmap")
#' plot(cluster_ex,type="chain")
#' 
#' # Comparison of cluster configurations #
#' # 1) Minimum distance with average MCMC iterations
#' # 2) Minimum Heterogeneity Measure (HM)
#' plot( x=jitter(cluster_ex$cluster),y=jitter(cluster_ex$clusterHMmin), col="#FF000080", pch=20,
#'       main=paste("Comparison of two relevant cluster configurations"),
#'       xlab="minimizes distance to average MCMC grouping", ylab="minimizes Heterogeneity Measure" )
#' 
#' # Comparison with the original clusters in the simulated data
#' plot(x=jitter(Z_latent_ex_5_1$cluster),
#'      y=jitter(cluster_ex$cluster),
#'      main=paste("Comparison real configuration with the model results"),
#'      xlab="Real cluster",
#'      ylab="Model cluster",
#'      pch=19, col="#FF000080")
#' }
#'
#' ##############################
#' #      Households data       #
#' #    Carmona et al. (2017)   #
#' ##############################
#' 
#' # Testing "MIXclustering" function with poverty.data #
#' # Data and parameters are discussed in section 5.3 of Carmona et al. (2017) #
#' 
#' # Set seed for reproducibility #
#' set.seed(0) 
#' 
#'
#' \dontrun{
#' # relevant variables for clustering households #
#' Y_names <- c("ict_norm",
#'              "ic_ali","ic_asalud","ic_cv",
#'              "ic_rezedu","ic_sbv","ic_segsoc",
#'              "niv_ed","tam_loc")
#' Y_var_type <- c("c","o","o","o","o","o","o","o","m")
#' 
#' # using only data from state 15 (Edomex) #
#' aux_subset <- rep(T,nrow(poverty.data))
#' aux_subset <- aux_subset & is.element(substr(poverty.data$folioviv,1,2),"15")
#'
#' Y_data <- poverty.data[aux_subset,Y_names]
#'
#' ### Sampling probability dependin on the scenario ###
#' # Scenario description in section 5.3 of Carmona et al. (2017) #
#' # Choose 1, 2 or 3 #
#' poverty_sampling_spec <- 3
#' 
#' 
#' if (poverty_sampling_spec == 1) {
#'   k <- 1
#'   sampling_prob_pov <- rep(1,nrow(Y_data))
#' } else if (poverty_sampling_spec == 2) {
#'   k <- 2 * mean(poverty.data[aux_subset,"factor_hog"])
#'   sampling_prob_pov <- 1/poverty.data[aux_subset,"factor_hog"]
#' } else if (poverty_sampling_spec == 3) {
#'   k <- 4 * mean(poverty.data[aux_subset,"factor_hog"])
#'   sampling_prob_pov <- 1/poverty.data[aux_subset,"factor_hog"]
#' }
#' 
#' cluster_poverty <- MIXclustering( Y=Y_data,
#'                                   var_type=Y_var_type,
#'                                   n_iter_out=1500,
#'                                   n_burn=200,
#'                                   n_thin=3,
#'                                   
#'                                   alpha = 0.5,
#'                                   d_0_a = 1, d_1_a = 1,
#'                                   d_0_b = 1, d_1_b = 1,
#'                                   
#'                                   eta = 2,
#'                                   kappa = 5,
#'                                   delta = 4,
#'                                   
#'                                   d_0_z = 2.1, d_1_z = 30,
#'                                   d_0_mu = 2.1, d_1_mu = 30,
#'                                   
#'                                   sampling_prob = k * sampling_prob_pov )
#' 
#' summary(cluster_poverty)
#' plot(cluster_poverty,type="heatmap")
#' plot(cluster_poverty,type="chain")
#' }
#'
#' @references
#' Carmona, C., Nieto-Barajas, L. E. & Canale, A. (2017). \emph{Model based approach for household clustering with mixed scale variables.} (preprint: \url{http://arxiv.org/abs/1612.00083})
#'
#' @importFrom stats aggregate as.formula cor rgamma rmultinom sd
#' @import plyr
#'
#' @useDynLib BNPMIXcluster
#' @importFrom Rcpp sourceCpp
#' @importFrom stats weighted.mean 
#' @export

MIXclustering <- function( Y,
                           var_type,
                           
                           n_iter_out=2000,
                           n_burn=100,
                           n_thin=2,
                           
                           a_fix=NULL,
                           alpha=0.5,
                           d_0_a=1, d_1_a=1,
                           
                           b_fix=NULL,
                           d_0_b=1, d_1_b=1,
                           eta=2,
                           
                           d_0_z=2.1, d_1_z=30,
                           kappa=5, delta=4,
                           
                           d_0_mu=2.1, d_1_mu=30,
                           
                           sampling_prob=NULL,
                           expansion_f=NULL,
                           
                           log_file=NULL,
                           keep_param_chains=FALSE ) {
  
  if( !any( c(is.matrix(Y),is.data.frame(Y)) ) ) { stop("Y has to be a Matrix or data frame, with each column representing a variable, and each row representing an individual") }
  
  USING_CPP = TRUE
  max.time=Inf
  
  #on.exit(browser())
  
  log_clock <- Sys.time()
  
  if( !is.null(log_file) ) {
    cat("\n
        ***** MIXclustering *****\n\n",
        "----- characteristics -----\n",
        "effective iterations = ",n_iter_out,"\n\n",
        "alpha=",alpha,"\n",
        "d_0_a=",d_0_a,"\n",
        "d_1_a=",d_1_a,"\n",
        "eta=",eta,"\n",
        "d_0_b=",d_0_b,"\n",
        "d_1_b=",d_1_b,"\n",
        "kappa=",kappa,"\n",
        "delta=",delta,"\n",
        "d_0_z=",d_0_z,"\n",
        "d_1_z=",d_1_z,"\n",
        "d_0_mu=",d_0_mu,"\n",
        "d_1_mu=",d_1_mu,"\n",
        "---------------------------\n",
        "Starting time:\n",as.character(log_clock),"\n\n",
        "---------------------------\n\n",
        "Iteration times:\n\n",
        "iter_i , elapsed_minutes , Sys.time\n",
        file=log_file)
  }
  
  dev_verbose = FALSE
  
  cl <- match.call()
  
  ######
  #   total number of iterations
  #####
  n_iter <- n_burn + 1 + (n_iter_out-1)*(n_thin)
  
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
  
  # Sorting Y columns
  if(!is.null(colnames(Y))) {
    Y_colnames_orig <- colnames(Y)
  } else {
    colnames(Y) <- Y_colnames_orig <- paste("Y",1:ncol(as.matrix(Y)),sep=".")
  }
  Y_orig <- Y
  var_type_orig <- var_type
  
  Y_new_order <- c( which(is.element(var_type,var_type_all[1])),
                    which(is.element(var_type,var_type_all[2])),
                    which(is.element(var_type,var_type_all[3])))
  var_type <- var_type[Y_new_order]
  Y <- Y[,Y_new_order,drop=F]
  
  ##### Missing data #####
  # Only consider rows with complete information, i.e. no missing data allowed
  aux_na_Y <- apply(!is.na(Y),1,all)
  
  Y <- Y[aux_na_Y,,drop=F]
  if(sum(!aux_na_Y)>0){
    cat(sum(!aux_na_Y),' rows with missing data were removed\n')
  }
  
  # vector for returning Z to Y incorporating back the NAs
  map_z_to_orig <-cumsum(aux_na_Y)
  map_z_to_orig[duplicated(map_z_to_orig)]<-NA
  
  n <- nrow(Y)
  p <- ncol(Y)
  
  # continuous latent variables
  n_c <- sum( is.element( var_type, var_type_all[1] ) )
  # ordinal latent variables (categorical, ordered)
  n_o <- sum( is.element( var_type, var_type_all[2] ) )
  # nominal variables (categorical, unordered)
  n_m <- sum( is.element( var_type, var_type_all[3] ) )
  # nominal latent variables (categorical, unordered)
  n_m_l = rep(NA,n_m)
  if(n_m>0) {
    for (i in 1:n_m ) {
      # i=1
      aux_i = which( is.element( var_type, var_type_all[3] ) )[i]
      # k-1 latents are needed for a categorical variable with k categories
      n_m_l[i] = length(unique(Y[,aux_i]))-1
    }
    rm(i,aux_i)
  }
  # total number of latent variables in Z
  n_q <- n_c+n_o+sum(n_m_l)
  
  
  ### Transforms Y to numerical matrix ###
  Y_numerical <- matrix(NA,nrow(Y),ncol(Y))
  if( n_c>0 ) {
    for(j in which(is.element( var_type, var_type_all[1] ))) {
      Y_numerical[,j] <- as.numeric(Y[,j])
    }
  }
  if( n_o>0 ) {
    for(j in which(is.element( var_type, var_type_all[2] ))) {
      Y_numerical[,j] <- match(Y[,j],sort(unique(Y[,j])))
    }
  }
  if( n_m>0 ) {
    for(j in which(is.element( var_type, var_type_all[3] ))) {
      Y_numerical[,j] <- match(Y[,j],sort(unique(Y[,j])))
    }
  }
  Y <- Y_numerical
  rm(Y_numerical,j)
  
  # changes the colnames of Y for simplicity and standarization
  
  colnames(Y) <- paste("var_",
                       c(rep(var_type_all[1],n_c),rep(var_type_all[2],n_o),rep(var_type_all[3],n_m)),"_",
                       formatC( unlist(mapply(seq,1,c(n_c,n_o,n_m),length.out=c(n_c,n_o,n_m))) , width = 2, flag = '0'),
                       sep="")
  
  if(is.null(sampling_prob)) {
    
    if(is.null(expansion_f)) {
      expansion_f<-rep(1,n)
    } else {
      if(length(expansion_f)!=n) {
        cat('\nError: There is an inconsistency between "expansion_f" and the number of rows in "Y"\n')
        stop('There is an inconsistency between "expansion_f" and the number of rows in "Y"')
      }
    }
    
    sampling_prob <- 1/expansion_f # Sample design probabilities
    
    
    
    if(length(sampling_prob)!=n) {
      cat('\nError: There is an inconsistency between "sampling_prob" and the number of rows in "Y"\n')
      stop('There is an inconsistency between "sampling_prob" and the number of rows in "Y"')
    }
  } else {
    if(is.null(expansion_f)) {
    expansion_f <- 1/sampling_prob
    } else {
      stop('Only one of "sampling_prob" or "expansion_f" should be specified')
    }
  }
  
  #sampling_prob <- sampling_prob/sum(sampling_prob) # should add up 1
  
  if(n_c>0){
    ### scale continuos variables to have mean 0 and sd 1 ###
    Y[,1:n_c] <- scale(Y[,1:n_c],center=T,scale=T)
  }
  
  #####     Simulating latent variables 'Z' from 'Y'     #####
  if(FALSE) {
    # DANGER: This is not ready yet! #
    Z <- get_latents_cpp( Y=as.matrix(Y),
                          var_type = match(var_type,var_type_all),
                          mu_Z = matrix(0,nrow=n,ncol=n_q),
                          sigma_Z = diag(1,nrow=n_q,ncol=n_q),
                          Z_ini = matrix(0,1,1) )
  } else {
    Z <- get_latents( Y=Y,
                      var_type=var_type,
                      USING_CPP=USING_CPP )$Z
  }
  
  ##### Start: Initializing the chain #####
  
  # sigma_Z: variances in Lambda
  sigma_Z_diag <- apply(Z,2,sd)
  Lambda <- diag(sigma_Z_diag,nrow=length(sigma_Z_diag),ncol=length(sigma_Z_diag))
  rm(sigma_Z_diag)
  
  # unit variances for "sigma_Z"
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
  
  # mu_Z <- Z
  # it is not necessary to define an object for mu_Z, since it is is fully characterized by
  # its unique values 'mu_star'
  # and its mapping 'mu_star_map'
  
  # matrix with unique values of mu_Z
  if(n_c>0) {
    mu_star <- Z[!duplicated(Z[,1]),,drop=F]
  } else {
    mu_star <- Z[!duplicated(Z),,drop=F]
  }
  
  # mapping each mu_Z to mu_star
  mu_star_map <- match(data.frame(t(Z)), data.frame(t(mu_star)))
  
  # how many repeated values of each mu_star
  mu_star_n_r <- as.numeric(table(mu_star_map))
  
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
  
  if(keep_param_chains) {
    Lambda_sim <- matrix(as.numeric(NA),nrow=n_iter,ncol=ncol(sigma_Z) )
    Omega_sim <- array(as.numeric(NA),dim=c(nrow(sigma_Z),ncol(sigma_Z),n_iter))
    a_sim <- as.numeric(NULL)
    b_sim <- as.numeric(NULL)
    
    #Lambda_accept <- matrix(as.numeric(NA),nrow=n_iter,ncol=ncol(sigma_Z) )
    #Omega_accept <- array(as.numeric(NA),dim=c(nrow(sigma_Z),ncol(sigma_Z),n_iter))
    #a_accept <- as.numeric(NULL)
    #b_accept <- as.numeric(NULL)
  } else {
    Lambda_sim <- NULL
    Omega_sim <- NULL
    a_sim <- NULL
    b_sim <- NULL
    
    #Lambda_accept <- NULL
    #Omega_accept <- NULL
    #a_accept <- NULL
    #b_accept <- NULL
  }
  
  
  #if(dev_verbose) {
  cat("*** Clustering estimation started ***\n")
  #cat('     Sys.time() = ',as.character(Sys.time()),sep='')
  #}
  time_sim <- Sys.time()
  
  lsa_iter <- c(ls(),"iter_i","lsa_iter","pb")
  
  for( iter_i in 1:n_iter) {
    
    if(iter_i==1) {
      # initializing progress bar
      pb <- plyr::create_progress_bar("time")
      pb$init(n_iter)
      
      on.exit(pb$term())
    }
    
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
      v_i <- solve( solve(sampling_prob[i] * sigma_Z) + solve(sigma_mu) )
      # u_i <- v_i %*% solve(sampling_prob[i] * sigma_Z) %*% matrix(Z[i,],nrow=n_q,ncol=1)
      u_i <- v_i %*% solve(sampling_prob[i] * sigma_Z) %*% matrix(mu_star[mu_star_map[i],],nrow=n_q,ncol=1)
      
      # eliminates the non-simmetry due to numerical precision
      if( max( abs(v_i-t(v_i)) )<1e-7 ) {v_i <- v_i - (v_i-t(v_i))/2} else {browser()}
      
      # number of diferent values of mu_star[-1]
      mu_star_n_r_temp <- mu_star_n_r
      mu_star_n_r_temp[mu_star_map[i]] <- mu_star_n_r_temp[mu_star_map[i]]-1
      r_i <- sum(mu_star_n_r_temp>0)
      
      D_0 <- ( b + a * r_i ) * dmvnrm_arma( x=Z[i,,drop=F] , mean=rep(0,n_q), sigma=sampling_prob[i]*sigma_Z+sigma_mu )
      D_values <- as.numeric(NULL)
      for(r in 1:length(mu_star_n_r_temp)) {
        D_values[r] <- ( mu_star_n_r_temp[r] - a ) * dmvnrm_arma( x=Z[i,,drop=F] , mean=mu_star[r,], sigma=sampling_prob[i]*sigma_Z )
      }
      D_values[D_values<0] <- 0
      
      if(sum( D_0, D_values )==0) {
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
      if(any(is.na(mu_star_new))){
        stop("There is a problem sampling from mu_Z (code 1)")
      }
    } #finished simulating from mu_Z
    
    # Identifies values of mu_star that have no Z elements pointing to them
    aux_no_map <- is.element(mu_star_n_r_new,0)
    
    # Replaces new simulated values of mu_star
    mu_star <- mu_star_new[!aux_no_map,,drop=F]
    mu_star_map <- mu_star_map_new - cumsum(aux_no_map)[mu_star_map_new]
    mu_star_n_r <- mu_star_n_r_new[!aux_no_map]
    
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
        v_i_star <- solve(sum(1/sampling_prob[I_j]) * sigma_Z_inv + sigma_mu_inv )
      } else {
        for( k in 1:sum(I_j) ) {
          if(k==1) {
            v_i_star <- (1/sampling_prob[I_j])[k] * sigma_Z_inv + sigma_mu_inv
          } else {
            v_i_star <- v_i_star + (1/sampling_prob[I_j])[k] * sigma_Z_inv + sigma_mu_inv
          }
        }
        v_i_star <- solve(v_i_star)
        rm(k)
      }
      weighted_z <- apply( matrix(1/sampling_prob[I_j],nrow=sum(I_j),ncol=n_q,byrow=F) * matrix(Z[I_j,],nrow=sum(I_j),ncol=n_q), 2, sum )
      weighted_z <- matrix(weighted_z,nrow=n_q,ncol=1)
      u_i_star <- v_i_star %*% solve(sigma_Z) %*% weighted_z
      
      mu_star_new[j,] <- mvtnorm::rmvnorm(1,mean=u_i_star,sigma=v_i_star)
    }
    if(dev_verbose) {
      cat('...Done! \n')
    }
    # Replaces new simulated values of mu_star
    mu_star <- mu_star_new
    
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
    
    ### (d) Sampling sigma_Z, variances in Lambda ###
    if(dev_verbose) {
      cat( 'Sampling sigma_Z: variances in Lambda \n' )
    }
    
    # all elements of sigma_Z with variance different of 1
    for( j_sigma in which(!aux_var1_Z)) {
      # j_sigma <- 1
      #for( j_sigma in 1:n_q) {
      
      if(dev_verbose) {
        cat(j_sigma,", ")
      }
      
      aux_Lambda <-  sampling_Lambda_jj( n_sim_mh=1, sigma_jj_ini=Lambda_new[j_sigma,j_sigma]^2,j=j_sigma,
                                         d_0_z=d_0_z, d_1_z=d_1_z, kappa=kappa,
                                         Z=Z, mu_Z=mu_star[mu_star_map,,drop=F], sigma_Z=sigma_Z_new, sampling_prob=sampling_prob,
                                         max.time=max.time,n.burn=0,
                                         verbose=F,
                                         USING_CPP=USING_CPP )
      
      Lambda_new[j_sigma,j_sigma] <- sqrt(aux_Lambda[[1]])
      #Lambda_accept[iter_i,j_sigma] <- aux_Lambda[[2]]
      
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
    }
    
    if(dev_verbose) {
      cat('...Done! \n')
    }
    
    
    ### (e) Sampling sigma_Z, correlations in Omega ###
    if ( all(dim(Omega)>c(1,1)) ) {
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
          
          aux_omega_ij_new <- sampling_Omega_ij( n=1,
                                                 Omega.ini=Omega_new,
                                                 i=i_omega,
                                                 j=j_omega,
                                                 delta=delta,
                                                 Z=Z,
                                                 mu_Z=mu_star[mu_star_map,,drop=F],
                                                 Lambda=Lambda_new,
                                                 sampling_prob=sampling_prob,
                                                 n.burn=0,
                                                 n.thin=0,
                                                 max.time=max.time,
                                                 USING_CPP=USING_CPP )
          
          
          omega_ij_new <- aux_omega_ij_new[[1]]
          #Omega_accept[i_omega,j_omega,iter_i] <- Omega_accept[j_omega,i_omega,iter_i] <- aux_omega_ij_new[[2]]
          
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
      
      # updates 'sigma_Z' after all element in 'Omega' are sampled
      # only once because 'Omega' sampling IS NOT dependent of 'sigma_Z'
      sigma_Z_new <- Lambda_new %*% Omega_new %*% Lambda_new
      
      if( max( abs(sigma_Z_new-t(sigma_Z_new)) ) < 1e-7 ) {
        sigma_Z_new <- sigma_Z_new - (sigma_Z_new-t(sigma_Z_new))/2
      } else {
        stop("There is a problem sampling from sigma_Z (code 1)")
      }
      
      if( !matrixcalc::is.positive.definite(sigma_Z_new) ) {
        cat("*****\nProcess finished because 'sigma_Z_new' is not positive definite!\n*****");
        return()
      }
      
      # eliminates the non-simmetry due to numerical precision
      if( max( abs(sigma_Z_new-t(sigma_Z_new)) ) < 1e-7 ) {
        sigma_Z_new <- sigma_Z_new - (sigma_Z_new-t(sigma_Z_new))/2
      } else {
        stop("There is a problem sampling from sigma_Z (code 2)")
      }
      
      #if( !matrixcalc::is.positive.definite(sigma_Z_new) ) {
      #  cat("*****\nProcess finished because 'sigma_Z' is not positive definite!\n*****");
      #  return()
      #}
      
      sigma_Z <- sigma_Z_new
      Lambda <- Lambda_new
      Omega <- Omega_new
      if(dev_verbose) {
        cat('...Done! \n')
      }
    }
    
    # (f) Sampling "a"
    
    if(is.null(a_fix)){
      
      if(dev_verbose) {
        cat( 'Sampling "a":\n' )
      }
      
      aux_a_new <- sampling_a( n=1, a.ini=a,
                               b=b, alpha=alpha, d_0_a=d_0_a, d_1_a=d_1_a,
                               mu_star_n_r=mu_star_n_r,
                               max.time=max.time, n.burn=0,
                               verbose=F,
                               USING_CPP=USING_CPP)
      
      a_new <- aux_a_new$a.chain
      #a_accept <- c( a_accept , aux_a_new$accept.indic )
      
    } else {
      a_new <- a_fix
      #a_accept <- NULL
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
                               max.time=max.time, n.burn=0, eta=eta,
                               verbose=F,
                               USING_CPP=USING_CPP )
      
      b_new <- aux_b_new[[1]]
      #b_accept <- c(b_accept, aux_b_new[[2]])
      
      if(!all(b_new>-a)){stop('There is a problem sampling from "b", it should be >-a\nb=',b,"\n-a=",-a,sep="")}
      
    } else {
      b_new <- b_fix
      #b_accept <- NULL
    }
    b <- b_new
    if(dev_verbose) {
      cat('...Done! \n')
    }
    
    # (h) Sampling "Z_ij"
    
    if(dev_verbose) {
      cat( 'Sampling "Z_ij":\n' )
    }
    
    if( FALSE ) {
      # DANGER: This is not ready yet! #
      Z_new <- get_latents_cpp( Y = as.matrix(Y),
                                var_type = match(var_type, var_type_all),
                                mu_Z = mu_star[mu_star_map,,drop=F],
                                sigma_Z = sigma_Z,
                                Z_ini = Z )
    } else {
      Z_new <- get_latents( Y=Y,
                            var_type=var_type,
                            mu_Z=mu_star[mu_star_map,,drop=F],
                            sigma_Z=sigma_Z,
                            Z_old=Z,
                            USING_CPP=USING_CPP)$Z
    }
    if(dev_verbose) {
      cat('...Done! \n')
    }
    
    ### Storing simulation values for each iteration ###
    
    #mu_Z_sim[,,iter_i] <- mu_star[mu_star_map,]
    #sigma_Z_sim[,,iter_i] <- sigma_Z
    mu_star_map_sim[,paste("iter_",iter_i,sep="")] <- mu_star_map
    #mu_star_n_r_sim[[iter_i]] <- mu_star_n_r
    
    if(keep_param_chains) {
      Lambda_sim[iter_i,] <- diag(Lambda)
      Omega_sim[,,iter_i] <- Omega
      a_sim[iter_i] <- a
      b_sim[iter_i] <- b
    }
    
    time_sim[iter_i+1] <- Sys.time()
    
    if(dev_verbose) {
      cat('\n   Finished iter_i=',iter_i,', \n     elapsed time:\n          ',format(diff(time_sim)[iter_i]),' this iteration\n          ',format(diff(time_sim[c(1,length(time_sim))])), ' total\n**----**',sep='')
    }
    
    rm(list=setdiff(ls(),lsa_iter))
    pb$step()
    
    if( !is.null(log_file) & ( ((iter_i%%100)==0) || (iter_i==n_iter) ) ) {
      cat( as.character(iter_i)," , ",difftime(time1=Sys.time(),time2=log_clock,units="mins")," , ",as.character(Sys.time())," \n",
           file=log_file, append=TRUE )
    }
    
  }
  
  #####
  #   Clustering results summary   #
  #####
  
  sort.cluster <- function(x) {
    # function for sorting clustering numbers
    # cluster 1 will be the one with more elements
    plyr::mapvalues(x,
                    from=as.numeric(names(sort(table(x),decreasing=T))),
                    to=sort(unique(x)))
  }
  
  ### iterations that will be reported ###
  # after burn-in period and thinning
  iter_out <- seq(from=n_burn+1,to=n_iter,by=n_thin)
  
  ### reported clusters (after ordering) ###
  clusters <- sapply( mu_star_map_sim[,iter_out], sort.cluster )
  clusters <- data.frame(clusters)
  
  ### Number of clusters in each iteration ###
  n.clusters <- apply(clusters,2,max)
  
  ### Average cluster matrix ###
  cluster_matrix_sum <- matrix(0,nrow=dim(clusters)[1],ncol=dim(clusters)[1])
  for(iter_i in 1:dim(clusters)[2]) {
    # iter_i <- 1
    cluster_matrix_i <- matrix(0,nrow=dim(clusters)[1],ncol=dim(clusters)[1])
    aux_cluster_num <- clusters[,iter_i]
    
    # creates the clustering matrix for "iter_i"
    for(cluster_j in 1:max(aux_cluster_num)) {
      # cluster_j <- 1
      cluster_matrix_i[aux_cluster_num==cluster_j,aux_cluster_num==cluster_j] <- 1
    }
    
    # add clustering matrix with "cluster_matrix_sum"
    cluster_matrix_sum <- cluster_matrix_sum + cluster_matrix_i
  }
  cluster.matrix.avg <- cluster_matrix_sum/dim(clusters)[2]
  
  # Distance of each cluster to the average matrix #
  cluster_dist <- rep(0,dim(clusters)[2])
  for(iter_i in 1:dim(clusters)[2]) {
    #cat(iter_i,",")
    # iter_i <- 10
    cluster_matrix_i <- matrix(0,nrow=dim(clusters)[1],ncol=dim(clusters)[1])
    aux_cluster_num <- clusters[,iter_i]
    
    # creates the clustering matrix for "iter_i"
    for(cluster_j in 1:max(aux_cluster_num)) {
      # cluster_j <- 1
      cluster_matrix_i[aux_cluster_num==cluster_j,aux_cluster_num==cluster_j] <- 1
    }
    
    # calculates distance to the average matrix
    cluster_dist[iter_i] <- sum((cluster_matrix_i-cluster.matrix.avg)^2)
  }
  
  # Which cluster is closer to the average similarity matrix #
  cluster.cta <- clusters[,which.min(cluster_dist)]
  
  # Summary of original variables by cluster
  Y <- data.frame(cluster=cluster.cta,Y_orig) # restablishes original Y
  Y.cluster.summary <- data.frame(cluster=sort(unique(cluster.cta)),
                                  n.elements=as.numeric(table(cluster.cta)))
  
  for(i in 2:ncol(Y)) {
    if( is.element(var_type_orig[i-1],var_type_all[1]) ) {
      # continuous variables
      Y_summary_aux <- aggregate(as.formula(paste(colnames(Y)[i],"~cluster",sep="")),
                                 data=Y,
                                 FUN=mean,
                                 na.rm=T)
      Y.cluster.summary <- merge( Y.cluster.summary,Y_summary_aux,by="cluster",all.x=T)
    } else if( is.element(var_type_orig[i-1],var_type_all[-1]) ) {
      Y_summary_aux <- tapply(X=rep(1,nrow(Y)),INDEX=list(Y$cluster,Y[,i]),FUN=sum,na.rm=T)
      Y_summary_aux[is.na(Y_summary_aux)] <- 0
      Y_summary_aux <- Y_summary_aux / apply(Y_summary_aux,1,sum,na.rm=T)
      colnames(Y_summary_aux) <- paste(colnames(Y)[i],'value',colnames(Y_summary_aux),sep="_")
      Y_summary_aux <- data.frame(cluster=sort(unique(Y$cluster)),Y_summary_aux)
      Y.cluster.summary <- merge( Y.cluster.summary,Y_summary_aux,by="cluster",all.x=T)
    } else { stop('There is an unknown type in "var_type_orig"') }
  }
  
  ### Cluster Heterogeneity Measure ###
  q<-ncol(Y)
  hm <- vector("numeric")
  for( k in 1:dim(clusters)[2] ) {
    # k <- 1
    cluster<-clusters[,k]
    nc<-max(cluster)
    nr<-table(cluster)[1:nc]
    
    hmm<-matrix(NA,nc,q)
    
    for (i in 1:nc){
      # i <- 1
      for (j in 1:q){
        # j <- 1
        xbar<-weighted.mean(Y[cluster==i,j],expansion_f[cluster==i])
        x2bar<-weighted.mean(Y[cluster==i,j]^2,expansion_f[cluster==i])
        hmm[i,j]<-x2bar-xbar^2
      }
    }
    hm <- c( hm, sum(apply(hmm,1,sum)*nr) )
  }
  
  
  dpclust <- structure(
    list(
      cluster=cluster.cta[map_z_to_orig],
      clusterHMmin = clusters[map_z_to_orig,which.min(hm)],
      
      Y.cluster.summary=Y.cluster.summary,
      
      Y.var_type=var_type_orig,
      Y.na=which(!aux_na_Y),
      Y.n=n,
      Y.p=p,
      
      MC.clusters=clusters[map_z_to_orig,],
      
      # Heterogeneity Measure (HM)
      cluster.HM=hm[which.min(cluster_dist)],
      clusterHMmin.HM=min(hm),
      MC.clusters.HM=hm,
      
      cluster.matrix.avg=cluster.matrix.avg[map_z_to_orig,map_z_to_orig],
      
      #mu_star_n_r_sim=mu_star_n_r_sim,
      #mu_Z_sim=mu_Z_sim,
      #sigma_Z_sim=sigma_Z_sim,
      #mu_star_probs=mu_star_probs,
      
      MC.values =list(
        n.clusters=n.clusters,
        a=a_sim,
        b=b_sim,
        Lambda=Lambda_sim,
        Omega=Omega_sim
      ) ,
      MC.accept.rate=NULL,
      call=cl
    )
    ,class="MIXcluster")
  
  return(dpclust)
  
}
