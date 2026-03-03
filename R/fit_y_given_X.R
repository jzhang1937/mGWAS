#' Fit y given X via lasso
#'
#' @param data A list containing a matrix X and vector y.
#'
#' @return A fitted E(y|X) function
#' @export
#'
fit_y_given_X_lasso <- function(data){
  # Extract X and y.
  X = data$X
  y = data$y
  
  # Number of predictors
  p <- ncol(X)
  # Fit cross-validated Lasso
  cv.glmnet.fit <- glmnet::cv.glmnet(X, y)
  
  # Get CV Lasso coefficients
  glmnet.coefs = stats::coef(cv.glmnet.fit, s = "lambda.min")
  Z <- glmnet.coefs[2:(p+1)]
  
  fit <- function(X) {
    glmnet.coefs[1]+ X %*% Z
  }
  return(fit)
  
}

#' Fit y given X via post-lasso
#'
#' @param data A list containing a matrix X and vector y.
#'
#' @return A fitted E(y|X) function
#' @export
#'
fit_y_given_X_post_lasso <- function(data){
  # Extract X and y.
  X = data$X
  y = data$y
  
  # Number of predictors
  p <- ncol(X)
  # Fit cross-validated Lasso
  cv.glmnet.fit <- glmnet::cv.glmnet(X, y)
  
  # Get CV Lasso nonzero predictors
  glmnet.coefs = stats::coef(cv.glmnet.fit, s = "lambda.min")
  # Z <- glmnet.coefs[2:(p+1)]
  nonzero = which(glmnet.coefs[2:(p+1)] != 0)
  
  # Fit OLS with only the nonzero predictors if non-degenerate, otherwise return
  # CV Lasso
  if (length(nonzero) > 0) {
    Xsub = X[,nonzero]
    post_lasso = stats::lm(y ~ Xsub)
    # Post-Lasso predictor
    # Check for degeneracy
    if (any(is.na(post_lasso$coefficients))) {
      Z <- glmnet.coefs[2:(p+1)]
      fit <- function(X) {
        glmnet.coefs[1]+ X %*% Z
      }
    } else {
      fit <- function(X) {
        post_lasso$coefficients[1]+
          X[,nonzero] %*% post_lasso$coefficients[2:(length(nonzero)+1)]
      }
    }
    
  } else {
    # Default to CV Lasso
    Z <- glmnet.coefs[2:(p+1)]
    fit <- function(X) {
      glmnet.coefs[1]+X %*% Z
    }
  }
  return(fit)
}

#' Fit y given X via lasso
#'
#' @param data A list containing a matrix X and vector y.
#'
#' @return A fitted E(y|X) function
#' @export
#'
fit_y_given_X_lasso_binomial <- function(data){
  # Extract X and y.
  X = data$X
  y = data$y
  
  # Number of predictors
  p <- ncol(X)
  # Fit cross-validated Lasso
  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, family = "binomial")
  
  # Get CV Lasso coefficients
  glmnet.coefs = stats::coef(cv.glmnet.fit, s = "lambda.min")
  Z <- glmnet.coefs[2:(p+1)]
  
  fit <- function(X) {
    natural <- glmnet.coefs[1]+ X %*% Z
    exp(natural) / (1 + natural)
  }
  return(fit)
  
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fit_y_given_X_ranger <- function(data) {
  y <- data$y
  X <- data$X
  df <- data.frame(y = y, X = X)
  
  # fit the model
  model <- ranger::ranger(y ~ ., data = df, num.threads = 2)
  
  fit <- function(X) {
    # Check if new_data is a data frame or matrix
    if (!is.data.frame(X) && !is.matrix(X)) {
      stop("new_data must be a data frame or matrix.")
    }
    
    # Ensure the new data has the same number of features
    # expected_features <- model$feature_names
    # if (!all(expected_features %in% colnames(new_data))) {
    #  stop("new_data must contain the same features as the training data.")
    #}
    
    # Reorder columns to match the training data
    #new_data <- new_data[, expected_features]
    
    # Convert new data to matrix if necessary
    new_data <- as.data.frame(X)
    colnames(new_data) <- model$forest$independent.variable.names
    
    # Make predictions
    preds <- predict(model, new_data)
    
    return(preds$predictions)
  }
  
  return(fit)
  
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fit_y_given_X_ranger_categorical <- function(data) {
  y <- data$y
  X <- data$X
  df <- data.frame(y = y, X = X)
  
  # max categorical (must be an integer)
  K <- as.integer(max(df$y))
  df$y <- factor(df$y, levels = 0:K)
  
  # fit the model
  model <- ranger::ranger(y ~ ., data = df, probability = TRUE, num.threads = 2)
  
  fit <- function(X) {
    # Check if new_data is a data frame or matrix
    if (!is.data.frame(X) && !is.matrix(X)) {
      stop("new_data must be a data frame or matrix.")
    }
    
    # Ensure the new data has the same number of features
    # expected_features <- model$feature_names
    # if (!all(expected_features %in% colnames(new_data))) {
    #  stop("new_data must contain the same features as the training data.")
    #}
    
    # Reorder columns to match the training data
    #new_data <- new_data[, expected_features]
    
    # Convert new data to matrix if necessary
    new_data <- as.data.frame(X)
    colnames(new_data) <- model$forest$independent.variable.names
    
    # Make predictions
    preds <- predict(model, new_data)
    
    return(preds$predictions)
  }
  
  return(fit)
  
}

#' Title
#'
#' @param data
#'
#' @return
#' @export
#'
#' @examples
fit_y_given_X_ranger_binary <- function(data) {
  # Binary classification setup
  y <- data$y
  X <- data$X
  df <- data.frame(y = y, X = X)
  
  # Ensure binary factor levels ("0", "1")
  df$y <- factor(y, levels = c(0, 1))
  
  # Fit ranger model with probability output
  model <- ranger::ranger(y ~ ., data = df, probability = TRUE, num.threads = 2)
  
  # Prediction function: returns probability of class "1"
  fit <- function(X) {
    if (!is.data.frame(X) && !is.matrix(X)) {
      stop("X must be a data frame or matrix.")
    }
    
    new_data <- as.data.frame(X)
    colnames(new_data) <- model$forest$independent.variable.names
    
    preds <- predict(model, new_data)
    prob_1 <- preds$predictions[, "1"]  # probability of class "1"
    
    return(prob_1)
  }
  
  return(fit)
  
}

