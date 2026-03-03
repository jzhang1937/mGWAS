############################################################################################
#' \code{spa_cdf} SPA to CDF of T_n = S_n / sqrt(n)
#'
#' @param X Numeric vector of length n for the predictor variable.
#' @param Y Numeric vector of length n for the response variable.
#' @param X_on_Z_fit_vals X_on_Z_fit$fitted.values
#' @param Y_on_Z_fit_vals Y_on_Z_fit$fitted.values
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#' @param R stats::uniroot() search space endpoint
#' should be broadened.
#'
#' @return \code{test statistic}, \code{left-sided p-value}, \code{right-sided p-value},
#' \code{both-sided p-value}, and \code{spa.success} which specifies whether the saddlepoint
#' equation could be solved. If not, a backup method (GCM) had to be employed as a backup.
#'
#' @examples
#' n <- 100; p <- 2; normalize <- FALSE; return_cdf <- FALSE
#' data <- list(X = rbinom(n = n, size = 1, prob = 0.2),
#'              Y = rpois(n = n, lambda = 1),
#'              Z = matrix(rnorm(n = n*p, mean = 0, sd = 1), nrow = n, ncol = p))
#' X <- data$X; Y <- data$Y; Z <- data$Z
#' X_on_Z_fit <- suppressWarnings(stats::glm(X ~ Z, family = "binomial"))
#' Y_on_Z_fit <- suppressWarnings(stats::glm(Y ~ Z, family = "poisson"))
#' spacrt:::spa_cdf(X = X, Y = Y,
#'                  X_on_Z_fit_vals = X_on_Z_fit$fitted.values,
#'                  Y_on_Z_fit_vals = Y_on_Z_fit$fitted.values,
#'                  fam = "binomial", R = 100)
#'
#' @keywords internal
spa_cdf <- function(X, Y,
                    X_on_Z_fit_vals,
                    Y_on_Z_fit_vals,
                    fam,
                    R = 5){
  
  P <- X_on_Z_fit_vals
  W <- Y - Y_on_Z_fit_vals
  n <- length(P)
  
  # compute the products of residuals for each observation
  prod_resids <- (X - P) * W
  
  # compute the test statistic
  test_stat <- 1/sqrt(n) * sum(prod_resids)
  
  t <- test_stat + 1/sqrt(n) * sum(P*W)
  
  success_uniroot <- FALSE
  
  tryCatch({
    # try to solve the saddlepoint equation
    s.hat <- stats::uniroot(
      f = function(s) d1_wcgf(s, P = P, W = W, fam) - sqrt(n)*t,
      lower = -abs(R), upper = abs(R),
      extendInt = "yes",
      tol = .Machine$double.eps)$root
    
    success_uniroot <- TRUE
  }, error = function(e) {message("stats::uniroot() failed: ", conditionMessage(e))}
  )
  
  
  if(success_uniroot == TRUE && {
    suppressWarnings({
      r.hat <- sign(s.hat) * sqrt(2 * (sqrt(n)* s.hat * t -
                                       wcgf(s = s.hat, P = P, W = W, fam)))
      
      # Lugannani-Rice formula
      p.left <- stats::pnorm(r.hat) + stats::dnorm(r.hat) *
        (1/r.hat - 1/(s.hat*sqrt(d2_wcgf(s = s.hat, P = P, W = W, fam))))
    })
    
    # decide if p.left is NA or beyond the range [0, 1] or not
    all(p.left >= 0, p.left <= 1, !is.na(p.left))
  }
  ){
    res <- list(test_stat = t - 1/sqrt(n) * sum(P*W),
                p.left = p.left,
                p.right = 1 - p.left,
                p.both = 2*min(c(p.left, 1 - p.left)),
                spa.success = TRUE)
  }else {
    test_stat <- sum(prod_resids)/(stats::sd(prod_resids) * sqrt(n-1))
    
    res <- list(test_stat = test_stat,
                p.left = stats::pnorm(test_stat, lower.tail = TRUE),
                p.right = stats::pnorm(test_stat, lower.tail = FALSE),
                p.both = 2*stats::pnorm(abs(test_stat), lower.tail = FALSE),
                spa.success = FALSE)
  }
  
  return(res)
}


############################################################################################
#' \code{wcgf} is a function computing the cumulant generating function (CGF) of
#' distributions, multiplied by a weight function, from the GLM family
#'
#' @param s The point where the CGF will be computed.
#' @param P A vector containing the parameter values of the family of distributions.
#' @param W A vector containing the weights.
#' @param fam The GLM family which includes the distribution whose CGF is being
#' evaluated (values can be \code{gaussian}, \code{binomial}, \code{poisson}, etc).
#'
#' @return CGF of the weighted distribution evaluated at \code{s}.
#'
#' @keywords internal
wcgf <- function(s, P, W, fam){
  
  if(fam == 'binomial') return(sum(log(exp(s*W)*P + 1 - P)))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P*(exp(s*W) - 1)))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}



############################################################################################
#' \code{d1_wcgf} is a function computing the derivative of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return First derivative of CGF of the weighted distribution evaluated at \code{s}.
#'
#' @keywords internal
d1_wcgf <- function(s, P, W, fam){
  
  # if(fam == 'binomial') return(sum((W*P*exp(s*W)) / (exp(s*W)*P + 1 - P)))
  if(fam == 'binomial') return(sum((W*P) / (P + (1 - P) * exp(-s*W))))
  if(fam == 'gaussian') return()
  if(fam == 'Gamma') return()
  if(fam == 'inverse.gaussian') return()
  if(fam == 'poisson') return(sum(P * exp(s*W) * W))
  if(fam == 'quasi') return()
  if(fam == 'quasibinomial') return()
  if(fam == 'quasipoisson') return()
}


############################################################################################
#' \code{d2_wcgf} is a function computing the hessian of the weighted cumulant
#' generating function (WCGF) of distributions from GLM family
#'
#' @inheritParams wcgf
#'
#' @return Second derivative of CGF of the weighted distribution evaluated at \code{s}.
#'
#' @keywords internal
d2_wcgf <- function(s, P, W, fam){
  
  if(fam == 'binomial'){
    Q <- 1 - P
    # return(sum((W^2*P*Q*exp(s*W)) / (exp(s*W)*P + Q)^2))
    return(sum((W^2*P*Q) / (exp(s*W)*P^2 + 2 * P * Q + Q^2 * exp(-s*W))))
  }
  
  if(fam == 'gaussian'){
    Q <- 1 - P
    return()
  }
  
  if(fam == 'Gamma'){
    Q <- 1 - P
    return()
  }
  
  if(fam == 'inverse.gaussian'){
    Q <- 1 - P
    return()
  }
  
  if(fam == 'poisson'){
    return(sum(P*exp(s*W)*W^2))
  }
  
  if(fam == 'quasi'){
    Q <- 1 - P
    return()
  }
  
  if(fam == 'quasibinomial'){
    Q <- 1 - P
    return()
  }
  
  if(fam == 'quasipoisson'){
    Q <- 1 - P
    return()
  }
}
