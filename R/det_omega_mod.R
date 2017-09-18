
# MH Sampling from 'omega' #
det_omega_mod <- function(X,x_new,i,j) {
  X[i,j] <- X[j,i] <- x_new
  return(det(X))
}
