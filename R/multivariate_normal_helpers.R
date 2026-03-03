### Function to perform inversion of a submatrix using Woodbury Formula ###
get_submatrix_inverse <- function(A_inverse, j) {
  n <- nrow(A_inverse)
  
  # Extract the relevant blocks
  A_jj_inv <- A_inverse[j, j]
  A_j_inv <- A_inverse[j, -j]
  A_inv_j <- A_inverse[-j, j]
  A_inv_jj <- A_inverse[-j, -j]
  
  # Calculate the submatrix inverse using Woodbury formula
  submatrix_inverse <- A_inv_jj - A_inv_j %*% solve(A_jj_inv) %*% A_j_inv
  
  return(submatrix_inverse)
}

### Helper function that generates a function that resamples from a multivariate Gaussian ###
# Sigma is the covariance matrix
# prec is the precision matrix
# j is the index of the predictor
# B is number of resamples
# X is the data matrix to generate resamples for
generate_resamples_gaussian = function(mu, Sigma, prec, seed = 1234) {
  force(mu)
  force(Sigma)
  force(prec)
  generate_resamples = function(j, Xnew, B) {
    l <- nrow(Xnew)
    submat_inv = get_submatrix_inverse(prec, j)
    # vector to compute conditional mean (1 x p-1) transposed.
    condMeanOp <- t(Sigma[j,-j]%*%submat_inv)
    
    # conditional variance
    condVar <- max(0,Sigma[j,j] - Sigma[j,-j]%*%submat_inv%*%Sigma[-j,j])
    
    # Resample the jth coordinate for all subjects B times.
    resamp <- stats::rnorm(l*B,mu[j]+
                             sweep(Xnew[,-j], 2, mu[-j], FUN = '-')%*%condMeanOp,sqrt(condVar)) |>
      R.utils::withSeed(seed = seed)
    
    # Reshape the vector into a matrix with l rows and B columns
    resamp <- matrix(resamp, nrow = l)
    return(resamp)
  }
  return(generate_resamples)
}

### Helper function that generates function that computes conditional mean ###
# Sigma is the covariance matrix
# prec is the precision matrix
# j is the index of the predictor
# X is the data matrix to get conditional means for
conditional_mean_gaussian = function(mu, Sigma, prec) {
  force(mu)
  force(Sigma)
  force(prec)
  
  conditional_mean = function(j,Xnew) {
    submat_inv = get_submatrix_inverse(prec, j)
    # vector to compute conditional mean (p-1 x 1)
    # condMeanOp <- t(Sigma[j,-j]%*%prec[-j,-j])
    condMeanOp <- Sigma[j,-j]%*%submat_inv
    # Calculate conditional mean for all units
    # condMean <- mu[j]+
    #  sweep(Xnew[,-j], 2, mu[-j], FUN = '-')%*%condMeanOp
    condMean <- mu[j]+ condMeanOp %*%
      t(sweep(Xnew[,-j], 2, mu[-j], FUN = '-'))
    # return conditional mean; length l vector
    return(condMean)
    
  }
  return(conditional_mean)
  
}