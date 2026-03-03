#' Fit L(X) using knockoffs default
#' Code adapted from knockoff::create.second_order()
#'
#' @param shrink shrink or not
#' @param data An object contaning a design matrix X.
#'
#' @return An object with fields \code{generate_resamples} and \code{conditional_mean}
#' @export
fit_X_shrinkage <- function (data, shrink = T)
{
  # Extract the design matrix.
  X = data$X
  
  # Force X as it is used to make the manufactured functions.
  force(X)
  
  mu = colMeans(X)
  if (!shrink) {
    Sigma = stats::cov(X)
    if (!is_posdef(Sigma)) {
      shrink = TRUE
    }
  }
  if (shrink) {
    if (!requireNamespace("corpcor", quietly = T))
      stop("corpcor is not installed", call. = F)
    Sigma = tryCatch({
      suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,
                                                             verbose = F)), nrow = ncol(X)))
    }, warning = function(w) {
    }, error = function(e) {
      stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
    }, finally = {
    })
  }
  # estimated precision matrix
  prec <- suppressWarnings(matrix(as.numeric(corpcor::invcov.shrink(X,
                                                                 verbose = F)), nrow = ncol(X)))
  
  # Compute functions for generating resamples and conditional mean.
  generate_resamples = generate_resamples_gaussian(mu, Sigma, prec)
  conditional_mean = conditional_mean_gaussian(mu, Sigma, prec)
  
  return(list(generate_resamples = generate_resamples, conditional_mean = conditional_mean,
              prec = prec, sig = Sigma))
}