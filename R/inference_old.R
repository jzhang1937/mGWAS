#' Title
#'
#' @param data a list containing a vector of outcomes y, data matrix X, and optionally a tree
#' @param L number of modeled causal effects
#' @param K a similarity matrix
#' @param phylo whether to use the given tree to construct the similarity matrix K
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param method one of {'moments','MLE'} (sigma^2,tau^2) are estimated using method-of-moments or MLE
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns List containing PIP -- p x L matrix of PIPs, individually for each effect mu -- p x L matrix of posterior means conditional on causal omega -- p x L matrix of posterior precisions conditional on causal lbf -- length-L array of log-Bayes-factor for each effect lbf_variable -- p x L matrix of per-variable log-Bayes-factors ssq -- length-L array of final effect size variances s^2 sigmasq -- final value of sigma^2 tausq -- final value of tau^2 alpha -- length-p array of posterior means of infinitesimal effects
#' @export
#'
#' @examples
susieK_naive <- function(data, L, K = NULL, phylo = TRUE,
                         est_ssq = TRUE, ssq = NULL, ssq_range = c(0,1), pi0 = NULL,
                         est_sigmasq = TRUE, est_tausq = TRUE,
                         sigmasq = 1, tausq = 0,
                         method = "moments",
                         sigmasq_range = NULL, tausq_range = NULL,
                         PIP = NULL, mu = NULL,
                         maxiter = 100, PIP_tol = 1e-3, verbose = FALSE) {
  y <- data$y
  y <- y - mean(y)
  X <- data$X
  X <- sweep(X, 2, colMeans(X), FUN = "-")
  X <- sweep(X, 2, apply(X, 2, function(x) sqrt(mean((x - mean(x))^2))), FUN = "/")
  tree <- data$tree
  n <- nrow(X); p <- ncol(X)
  
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) stop("Sample names in y must match rownames of X")
  
  samples  <- names(y)
  variants <- colnames(X)
  
  if (is.null(K) & !phylo) {
    K <- X %*% t(X)
    message("K required if not using tree; using XX^T")
  }
  if (phylo) K <- ape::vcv.phylo(tree)[samples, samples]
  
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, p, L)
  if (is.null(mu))  mu  <- matrix(0,   p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf <- numeric(L)
  
  logpi0 <- if (is.null(pi0)) rep(log(1/p), p) else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  eig     <- eigen(K, symmetric = TRUE)
  eigvals <- rev(eig$values)
  Q       <- eig$vectors[, rev(seq_len(ncol(eig$vectors)))]
  # eigvals <- eig$values
  # Q       <- eig$vectors
  
  Z       <- crossprod(Q, X)   # Q^T X  (n x p)
  y_tilde <- crossprod(Q, y)   # Q^T y
  
  w <- 1 / (tausq * eigvals + sigmasq)
  
  diagXtOmegaX <- colSums(Z^2 * w)
  # diagXtOmegaX <- as.numeric(crossprod(Z^2, w))
  # diagXtOmegaX <- colSums(sweep(Z^2, 1, w, `*`))
  XtOmegay     <- crossprod(Z, w * y_tilde)
  
  # --- precompute invariants for MoMK ---
  yty  <- as.numeric(crossprod(y))
  Xty  <- crossprod(X, y)            # p-vector
  XtX  <- crossprod(X)               # p x p  (= Z^T Z since Q is orthogonal)
  XtXsq <- XtX %*% XtX              # p x p
  
  trXtX  <- sum(diag(XtX))
  # trXtKX <- sum(colSums(Z^2) * eigvals) # tr(X^T K X) = tr(Z^T diag(eigvals) Z)
  trXtKX <- sum(rowSums(Z^2) * eigvals)
  
  trK <- sum(eigvals)
  determinant_MoM <- det(matrix(c(n, trK, trXtX, trXtKX), 2, 2, byrow = TRUE))
  
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  
  for (it in seq_len(maxiter)) {
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    for (l in seq_len(L)) {
      b <- rowSums(mu * PIP) - mu[,l] * PIP[,l]
      
      XtOmegaXb <- crossprod(Z, w * (Z %*% b))
      XtOmegar  <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0
          -logsumexp(val)
        }
        opt   <- optimize(f, ssq_range, tol = 1e-5)
        ssq[l] <- opt$minimum
        if (verbose) cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
      }
      
      omega[,l]        <- diagXtOmegaX + 1 / ssq[l]
      mu[,l]           <- XtOmegar / omega[,l]
      lbf_variable[,l] <- XtOmegar^2 / (2 * omega[,l]) - 0.5 * log(omega[,l] * ssq[l])
      
      logPIP  <- lbf_variable[,l] + logpi0
      lbf[l]  <- logsumexp(logPIP)
      PIP[,l] <- exp(logPIP - lbf[l])
    }
    
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        res <- MoMK_naive(PIP, mu, omega, sigmasq, tausq, n,
                          XtX, XtXsq, Xty, yty,          # <-- updated arguments
                          trK = trK, trXtX, trXtKX,
                          est_sigmasq, est_tausq, verbose)
        sigmasq <- res$sigmasq
        tausq   <- res$tausq
      } else {
        stop("MLE not implemented for K")
      }
      
      w            <- 1 / (tausq * eigvals + sigmasq)
      diagXtOmegaX <- colSums(Z^2 * w)
      # diagXtOmegaX <- as.numeric(crossprod(Z^2, w))
      # diagXtOmegaX <- colSums(sweep(Z^2, 1, w, `*`))
      XtOmegay     <- crossprod(Z, w * y_tilde)
    }
    
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("Max PIP diff:", PIP_diff, "\n")
    if (PIP_diff < PIP_tol) { converged <- TRUE; break }
  }
  
  if (!exists("converged")) converged <- FALSE
  
  b          <- rowSums(mu * PIP)
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(PIP = PIP, marginalPIP = marginalPIP, mu = mu, omega = omega,
       lbf = lbf, lbf_variable = lbf_variable,
       ssq = ssq, sigmasq = sigmasq,
       tausq = tausq, converged = converged, variants = variants, determinant_MoM = determinant_MoM)
}

#' Title
#'
#' @param data a list containing a vector of outcomes y, data matrix X, and optionally a tree
#' @param L number of modeled causal effects
#' @param K a similarity matrix
#' @param phylo whether to use the given tree to construct the similarity matrix K
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param method one of {'moments','MLE'} (sigma^2,tau^2) are estimated using method-of-moments or MLE
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns List containing PIP -- p x L matrix of PIPs, individually for each effect mu -- p x L matrix of posterior means conditional on causal omega -- p x L matrix of posterior precisions conditional on causal lbf -- length-L array of log-Bayes-factor for each effect lbf_variable -- p x L matrix of per-variable log-Bayes-factors ssq -- length-L array of final effect size variances s^2 sigmasq -- final value of sigma^2 tausq -- final value of tau^2 alpha -- length-p array of posterior means of infinitesimal effects
#' @export
#'
#' @examples
susieK_old <- function(data, L, K = NULL, phylo = TRUE,
                   est_ssq = TRUE, ssq = NULL, ssq_range = c(0,1), pi0 = NULL,
                   est_sigmasq = TRUE, est_tausq = TRUE,
                   sigmasq = 1, tausq = 0,
                   method = "moments",
                   sigmasq_range = NULL, tausq_range = NULL,
                   PIP = NULL, mu = NULL,
                   maxiter = 100, PIP_tol = 1e-3, LD = FALSE, verbose = FALSE) {
  y <- data$y
  y <- y - mean(y)
  X <- data$X
  X <- sweep(X, 2, colMeans(X), FUN = "-")
  X <- sweep(X, 2, apply(X, 2, function(x) sqrt(mean((x - mean(x))^2))), FUN = "/")
  tree <- data$tree
  n <- nrow(X); p <- ncol(X)
  if (LD) {
    ld <- t(X) %*% X / n
  } else {
    ld = NULL
  }
  
  
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) stop("Sample names in y must match rownames of X")
  
  samples  <- names(y)
  variants <- colnames(X)
  
  if (is.null(K) & !phylo) {
    K <- X %*% t(X)
    message("K required if not using tree; using XX^T")
  }
  if (phylo) K <- ape::vcv.phylo(tree)[samples, samples]
  
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, p, L)
  if (is.null(mu))  mu  <- matrix(0,   p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf <- numeric(L)
  
  logpi0 <- if (is.null(pi0)) rep(log(1/p), p) else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  # --- eigendecomposition of K (n x n) for the GLS weights ---
  eig     <- eigen(K, symmetric = TRUE)
  eigvals <- rev(eig$values)
  Q       <- eig$vectors[, rev(seq_len(ncol(eig$vectors)))]
  
  Z       <- crossprod(Q, X)     # Q'X  (n x p)
  y_tilde <- crossprod(Q, y)     # Q'y  (n-vector)
  
  w <- 1 / (tausq * eigvals + sigmasq)
  
  diagXtOmegaX <- colSums(Z^2 * w)
  XtOmegay     <- crossprod(Z, w * y_tilde)
  
  # --- eigendecomposition of XX' (n x n) for fast MoMK ---
  eig_XX    <- eigen(tcrossprod(X), symmetric = TRUE)   # one-time O(n^3)
  D_Z       <- rev(eig_XX$values)
  U         <- eig_XX$vectors[, rev(seq_len(n))]
  Z_new     <- crossprod(U, X)                           # U'X  (n x p)
  y_new     <- as.numeric(crossprod(U, y))               # U'y  (n-vector)
  
  diagXtX   <- colSums(Z_new^2)                          # diag(X'X),   O(np)
  diagXtXsq <- as.numeric(crossprod(Z_new^2, D_Z))       # diag((X'X)^2), O(np)
  trXtX     <- sum(D_Z)                                  # tr(X'X) = sum of eigenvalues of XX'
  
  # --- other MoMK invariants ---
  yty  <- as.numeric(crossprod(y))
  Xty  <- crossprod(X, y)    # p-vector
  
  trK <- sum(eigvals)
  trXtKX <- sum(rowSums(Z^2) * eigvals)
  determinant_MoM <- det(matrix(c(n, trK, trXtX, trXtKX), 2, 2, byrow = TRUE))
  
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  
  converged <- FALSE
  
  for (it in seq_len(maxiter)) {
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    for (l in seq_len(L)) {
      b <- rowSums(mu * PIP) - mu[,l] * PIP[,l]
      
      XtOmegaXb <- crossprod(Z, w * (Z %*% b))
      XtOmegar  <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0
          -logsumexp(val)
        }
        opt    <- optimize(f, ssq_range, tol = 1e-5)
        ssq[l] <- opt$minimum
        if (verbose) cat(sprintf("Update s^2 for effect %d to %f\n", l, ssq[l]))
      }
      
      omega[,l]        <- diagXtOmegaX + 1 / ssq[l]
      mu[,l]           <- XtOmegar / omega[,l]
      lbf_variable[,l] <- XtOmegar^2 / (2 * omega[,l]) - 0.5 * log(omega[,l] * ssq[l])
      
      logPIP  <- lbf_variable[,l] + logpi0
      lbf[l]  <- logsumexp(logPIP)
      PIP[,l] <- exp(logPIP - lbf[l])
    }
    
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        res <- MoMK(PIP, mu, omega, sigmasq, tausq, n,
                    Z_new, D_Z, y_new,
                    Xty, yty,
                    diagXtX, diagXtXsq,
                    trK = trK, trXtX, trXtKX,
                    est_sigmasq, est_tausq, verbose)
        sigmasq <- res$sigmasq
        tausq   <- res$tausq
      } else {
        stop("MLE not implemented for K")
      }
      
      w            <- 1 / (tausq * eigvals + sigmasq)
      diagXtOmegaX <- colSums(Z^2 * w)
      XtOmegay     <- crossprod(Z, w * y_tilde)
    }
    
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("Max PIP diff:", PIP_diff, "\n")
    if (PIP_diff < PIP_tol) { converged <- TRUE; break }
  }
  
  n_iterations <- if (exists("it")) it else 0
  
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(PIP = PIP, marginalPIP = marginalPIP, mu = mu, omega = omega,
       lbf = lbf, lbf_variable = lbf_variable,
       ssq = ssq, sigmasq = sigmasq,
       tausq = tausq, converged = converged, 
       n_iterations = n_iterations, variants = variants,
       determinant_MoM = determinant_MoM, LD = ld)
}

# ──────────────────────────────────────────────────────────────────────────────
# susieKv2
#
# SuSiE with an arbitrary similarity matrix K for the residual covariance.
# Model:  Y = X beta + epsilon,   epsilon ~ N(0, tau^2 K + sigma^2 I)
#
# Identical interface to susieK; only the MoM step differs.
# ──────────────────────────────────────────────────────────────────────────────
#' Title
#' @param data a list containing a vector of outcomes y, data matrix X, and optionally a tree
#' @param L number of modeled causal effects
#' @param K a similarity matrix
#' @param phylo whether to use the given tree to construct the similarity matrix K
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns
#' @export
#'
#' @examples
susieKv2_old <- function(data, L, K = NULL, phylo = TRUE,
                     est_ssq     = TRUE,  ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                     est_sigmasq = TRUE,  est_tausq = TRUE,
                     sigmasq = 1, tausq = 0,
                     sigmasq_range = NULL, tausq_range = NULL,
                     PIP = NULL, mu = NULL,
                     maxiter = 100, PIP_tol = 1e-3, LD = FALSE, verbose = FALSE) {
  
  # ── data prep ────────────────────────────────────────────────────────────────
  y    <- data$y;   y <- y - mean(y)
  X    <- data$X
  X    <- sweep(X, 2, colMeans(X), FUN = "-")
  X    <- sweep(X, 2, apply(X, 2, function(x) sqrt(mean(x^2))), FUN = "/")
  tree <- data$tree
  
  n <- nrow(X); p <- ncol(X)
  if (LD) {
    ld <- t(X) %*% X / n
  } else {
    ld = NULL
  }
  
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_",  seq_len(n))
  if (is.null(names(y)))    names(y)    <- rownames(X)
  if (!all(names(y) %in% rownames(X))) stop("Sample names in y must match rownames of X")
  
  samples  <- names(y)
  variants <- colnames(X)
  
  if (is.null(K) & !phylo) {
    K <- X %*% t(X)
    message("K required if not using tree; using XX^T")
  }
  if (phylo) K <- ape::vcv.phylo(tree)[samples, samples]
  
  # ── initialisations ──────────────────────────────────────────────────────────
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1 / p, p, L)
  if (is.null(mu))  mu  <- matrix(0,     p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf          <- numeric(L)
  
  logpi0 <- if (is.null(pi0)) {
    rep(log(1 / p), p)
  } else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  # ── eigen decomposition of K ─────────────────────────────────────────────────
  eig     <- eigen(K, symmetric = TRUE)
  eigvals <- rev(eig$values)
  Q       <- eig$vectors[, rev(seq_len(ncol(eig$vectors)))]
  
  Z       <- crossprod(Q, X)   # Q'X,  n x p
  y_tilde <- crossprod(Q, y)   # Q'y,  length-n
  
  # ── precompute invariants ─────────────────────────────────────────────────────
  trK  <- sum(eigvals)
  trK2 <- sum(eigvals^2)
  
  yty  <- as.numeric(crossprod(y_tilde))           # y'y = ||y_tilde||^2
  ytKy <- sum(eigvals * y_tilde^2)                 # y'Ky
  
  Xty   <- as.numeric(crossprod(Z, y_tilde))       # X'y  = Z'y_tilde
  XtKy  <- as.numeric(crossprod(Z, eigvals * y_tilde))  # X'Ky = Z'(Lambda y_tilde)
  
  diagXtX  <- colSums(Z^2)                         # diag(X'X),  length-p
  diagXtKX <- as.numeric(crossprod(Z^2, eigvals))  # diag(X'KX), length-p
  
  # ── GLS weights and sufficient stats ─────────────────────────────────────────
  w            <- 1 / (tausq * eigvals + sigmasq)  # Omega^{-1} diagonal in eigenbasis
  diagXtOmegaX <- as.numeric(crossprod(Z^2, w))    # diag(X' Omega^{-1} X)
  XtOmegay     <- as.numeric(crossprod(Z, w * y_tilde))
  
  determinant_MoM <- det(matrix(c(n, trK, trK, trK2), 2, 2, byrow = TRUE))
  
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  
  # ── main loop ─────────────────────────────────────────────────────────────────
  converged <- FALSE
  
  for (it in seq_len(maxiter)) {
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    # -- single-effect updates --
    for (l in seq_len(L)) {
      b_minus_l     <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
      XtOmegaXb     <- as.numeric(crossprod(Z, w * as.numeric(Z %*% b_minus_l)))
      XtOmegar      <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0
          -logsumexp(val)
        }
        opt    <- optimize(f, ssq_range, tol = 1e-5)
        ssq[l] <- opt$minimum
        if (verbose) cat(sprintf("  s^2[%d] -> %f\n", l, ssq[l]))
      }
      
      omega[, l]        <- diagXtOmegaX + 1 / ssq[l]
      mu[, l]           <- XtOmegar / omega[, l]
      lbf_variable[, l] <- XtOmegar^2 / (2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
      
      logPIP  <- lbf_variable[, l] + logpi0
      lbf[l]  <- logsumexp(logPIP)
      PIP[, l] <- exp(logPIP - lbf[l])
    }
    
    # -- variance component update (MoMKv2) --
    if (est_sigmasq || est_tausq) {
      res <- MoMKv2(
        PIP = PIP, mu = mu, omega = omega,
        sigmasq = sigmasq, tausq = tausq, n = n,
        Z = Z, eigvals = eigvals, y_tilde = y_tilde,
        Xty = Xty, XtKy = XtKy,
        yty = yty, ytKy = ytKy,
        diagXtX = diagXtX, diagXtKX = diagXtKX,
        trK = trK, trK2 = trK2,
        est_sigmasq = est_sigmasq, est_tausq = est_tausq,
        verbose = verbose
      )
      sigmasq <- res$sigmasq
      tausq   <- res$tausq
      
      w            <- 1 / (tausq * eigvals + sigmasq)
      diagXtOmegaX <- as.numeric(crossprod(Z^2, w))
      XtOmegay     <- as.numeric(crossprod(Z, w * y_tilde))
    }
    
    # -- convergence check --
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("  Max PIP diff:", PIP_diff, "\n")
    if (PIP_diff < PIP_tol) { converged <- TRUE; break }
  }
  
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(
    PIP         = PIP,
    marginalPIP = marginalPIP,
    mu          = mu,
    omega       = omega,
    lbf         = lbf,
    lbf_variable = lbf_variable,
    ssq         = ssq,
    sigmasq     = sigmasq,
    tausq       = tausq,
    converged   = converged,
    variants    = variants,
    determinant_MoM = determinant_MoM,
    LD = ld
  )
}

# ──────────────────────────────────────────────────────────────────────────────
# susieKv3
#
# SuSiE with an arbitrary similarity matrix K.
# Model:  Y = X beta + epsilon,   epsilon ~ N(0, tau^2 K + sigma^2 I)
#
# Uses VCupdateK (ELBO-based IRLS) instead of MoM for (sigma^2, tau^2).
# The precomputation is also simpler: only the eigenbasis quantities are needed.
# ──────────────────────────────────────────────────────────────────────────────
susieKv3_old <- function(data, L, K = NULL, phylo = TRUE,
                     est_ssq = TRUE, ssq = NULL, ssq_range = c(0, 1), pi0 = NULL,
                     est_sigmasq = TRUE, est_tausq = TRUE,
                     sigmasq = 1, tausq = 0,
                     sigmasq_range = NULL, tausq_range = NULL,
                     PIP = NULL, mu = NULL,
                     maxiter = 100, PIP_tol = 1e-3,
                     maxiter_vc = 20, tol_vc = 1e-8,
                     LD = FALSE, 
                     verbose = FALSE) {
  
  # ── data prep ────────────────────────────────────────────────────────────────
  y    <- data$y;  y <- y - mean(y)
  X    <- data$X
  X    <- sweep(X, 2, colMeans(X), FUN = "-")
  X    <- sweep(X, 2, apply(X, 2, function(x) sqrt(mean(x^2))), FUN = "/")
  tree <- data$tree
  
  n <- nrow(X); p <- ncol(X)
  if (LD) {
    ld <- t(X) %*% X / n
  } else {
    ld = NULL
  }
  
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq_len(n))
  if (is.null(names(y)))    names(y)    <- rownames(X)
  if (!all(names(y) %in% rownames(X))) stop("Sample names in y must match rownames of X")
  
  samples  <- names(y)
  variants <- colnames(X)
  
  if (is.null(K) & !phylo) { K <- X %*% t(X); message("K required if not using tree; using XX^T") }
  if (phylo) K <- ape::vcv.phylo(tree)[samples, samples]
  
  # ── initialisations ───────────────────────────────────────────────────────────
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1 / p, p, L)
  if (is.null(mu))  mu  <- matrix(0,     p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf          <- numeric(L)
  
  logpi0 <- if (is.null(pi0)) rep(log(1 / p), p) else {
    out <- rep(-Inf, p); out[pi0 > 0] <- log(pi0[pi0 > 0]); out
  }
  
  # ── eigendecomposition of K ───────────────────────────────────────────────────
  eig     <- eigen(K, symmetric = TRUE)
  eigvals <- rev(eig$values)
  Q       <- eig$vectors[, rev(seq_len(ncol(eig$vectors)))]
  
  Z       <- crossprod(Q, X)            # Q'X,  n x p
  y_tilde <- as.numeric(crossprod(Q, y)) # Q'y, length-n
  
  # ── GLS weights and sufficient statistics ────────────────────────────────────
  w            <- 1 / (tausq * eigvals + sigmasq)
  diagXtOmegaX <- as.numeric(crossprod(Z^2, w))
  XtOmegay     <- as.numeric(crossprod(Z, w * y_tilde))
  
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  
  # ── main loop ─────────────────────────────────────────────────────────────────
  converged <- FALSE
  
  for (it in seq_len(maxiter)) {
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    # -- single-effect updates --
    for (l in seq_len(L)) {
      b_minus_l <- rowSums(mu * PIP) - mu[, l] * PIP[, l]
      XtOmegaXb <- as.numeric(crossprod(Z, w * as.numeric(Z %*% b_minus_l)))
      XtOmegar  <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0
          -logsumexp(val)
        }
        ssq[l] <- optimize(f, ssq_range, tol = 1e-5)$minimum
        if (verbose) cat(sprintf("  s^2[%d] -> %f\n", l, ssq[l]))
      }
      
      omega[, l]        <- diagXtOmegaX + 1 / ssq[l]
      mu[, l]           <- XtOmegar / omega[, l]
      lbf_variable[, l] <- XtOmegar^2 / (2 * omega[, l]) - 0.5 * log(omega[, l] * ssq[l])
      
      logPIP   <- lbf_variable[, l] + logpi0
      lbf[l]   <- logsumexp(logPIP)
      PIP[, l] <- exp(logPIP - lbf[l])
    }
    
    # -- variance component update (ELBO-based IRLS) --
    if (est_sigmasq || est_tausq) {
      res <- VCupdateK(PIP, mu, omega, sigmasq, tausq,
                       Z, eigvals, y_tilde,
                       est_sigmasq, est_tausq,
                       maxiter_vc, tol_vc, verbose)
      sigmasq <- res$sigmasq
      tausq   <- res$tausq
      
      w            <- 1 / (tausq * eigvals + sigmasq)
      diagXtOmegaX <- as.numeric(crossprod(Z^2, w))
      XtOmegay     <- as.numeric(crossprod(Z, w * y_tilde))
    }
    
    # -- convergence check --
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("  Max PIP diff:", PIP_diff, "\n")
    if (PIP_diff < PIP_tol) { converged <- TRUE; break }
  }
  
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(
    PIP          = PIP,
    marginalPIP  = marginalPIP,
    mu           = mu,
    omega        = omega,
    lbf          = lbf,
    lbf_variable = lbf_variable,
    ssq          = ssq,
    sigmasq      = sigmasq,
    tausq        = tausq,
    converged    = converged,
    variants     = variants,
    LD = ld
  )
}