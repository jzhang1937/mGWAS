jaccard_sim <- function(mat) {
  mat <- mat > 0
  inter <- tcrossprod(mat)  
  sums <- rowSums(mat)
  union <- outer(sums, sums, "+") - inter
  inter / union  
}

logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

MoM <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, Xty, yty,
                est_sigmasq, est_tausq, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  A <- matrix(c(n, sum(Dsq), sum(Dsq), sum(Dsq^2)), 2, 2)
  
  b <- rowSums(mu * PIP)
  # Vtb <- crossprod(V, b)
  Vtb  <- as.vector(crossprod(V, b))
  diagVtMV <- Vtb^2
  tmpD <- numeric(p)
  
  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    # Vtbl <- crossprod(V, bl)
    Vtbl <- as.vector(crossprod(V, bl))
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }
  
  diagVtMV <- diagVtMV + colSums(t(V^2) * tmpD)
  # diagVtMV <- diagVtMV + colSums(V^2 * tmpD)
  
  x0 <- yty - 2 * sum(b * Xty) + sum(Dsq * diagVtMV)
  x1 <- sum(Xty^2) - 2 * sum(Vtb * VtXty * Dsq) + sum(Dsq^2 * diagVtMV)
  x <- c(x0, x1)
  
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]; tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n; tausq <- 0
    }
    if (verbose) cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1,2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

MLE <- function(PIP, mu, omega, sigmasq, tausq, n, V, Dsq, VtXty, yty,
                est_sigmasq, est_tausq, sigmasq_range, tausq_range, it, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  if (is.null(sigmasq_range)) sigmasq_range <- c(0.2 * yty / n, 1.2 * yty / n)
  if (is.null(tausq_range)) tausq_range <- c(1e-12, 1.2 * yty / (n * p))
  
  b <- rowSums(mu * PIP)
  Vtb <- crossprod(V, b)
  
  diagVtMV <- Vtb^2
  tmpD <- numeric(p)
  
  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    Vtbl <- crossprod(V, bl)
    diagVtMV <- diagVtMV - Vtbl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }
  
  diagVtMV <- diagVtMV + colSums((t(V)^2) * tmpD)
  
  f <- function(x) {
    s2 <- x[1]; t2 <- x[2]
    denom <- t2 * Dsq + s2
    0.5 * (n - p) * log(s2) + 0.5 / s2 * yty +
      sum(0.5 * log(denom)
          - 0.5 * t2 / s2 * VtXty^2 / denom
          - Vtb * VtXty / denom
          + 0.5 * Dsq / denom * diagVtMV)
  }
  
  if (est_tausq) {
    opt <- optim(c(sigmasq, tausq), f, method = "L-BFGS-B",
                 lower = c(sigmasq_range[1], tausq_range[1]),
                 upper = c(sigmasq_range[2], tausq_range[2]))
    if (opt$convergence == 0) {
      sigmasq <- opt$par[1]; tausq <- opt$par[2]
      if (verbose) cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
    }
  } else if (est_sigmasq) {
    g <- function(x) f(c(x, tausq))
    opt <- optimize(g, sigmasq_range)
    sigmasq <- opt$minimum
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

MoMK_old <- function(PIP, mu, omega, sigmasq, tausq, n, eigvals,
                 XtX_diag, XTKX_diag, Xty, XTKy, yty, yKy, trK, trK2,
                 est_sigmasq, est_tausq, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  A <- matrix(c(n, trK, trK, trK2), 2, 2)
  
  b <- rowSums(mu * PIP)
  
  diagM <- b^2
  tmpD <- numeric(p)
  
  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    diagM <- diagM - bl^2
    tmpD <- tmpD + PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
  }
  
  diagM <- diagM + tmpD
  
  x0 <- yty - 2 * sum(b * Xty) + sum(XtX_diag * diagM)
  x1 <- yKy - 2 * sum(b * XTKy) + sum(XTKX_diag * diagM)
  
  x <- c(x0, x1)
  
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]; tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n; tausq <- 0
    }
    if (verbose) cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1,2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

MoMK_naive <- function(PIP, mu, omega, sigmasq, tausq, n,
                 XtX, XtXsq, Xty, yty,
                 trK, trXtX, trXtKX,
                 est_sigmasq, est_tausq, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  # A matches Python: [[n, trK], [trXtX, trXtKX]]
  A <- matrix(c(n, trK,
                trXtX, trXtKX), 2, 2, byrow = TRUE)
  
  b <- rowSums(mu * PIP)       # length-p posterior mean vector
  XtXb  <- XtX  %*% b
  XtXsqb <- XtXsq %*% b
  
  diag_XtX   <- diag(XtX)
  diag_XtXsq <- diag(XtXsq)
  
  # tr(X'X M) and tr(X'X^2 M)
  # = b'(X'X)b + sum_l [ diag(X'X)·s_l - bl'(X'X)bl ]
  trXtXM   <- as.numeric(b %*% XtXb)
  trXtXsqM <- as.numeric(b %*% XtXsqb)
  
  for (l in seq_len(L)) {
    bl <- mu[, l] * PIP[, l]
    sl <- PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
    XtXbl   <- XtX   %*% bl
    XtXsqbl <- XtXsq %*% bl
    # trXtXM   <- trXtXM   + sum(diag_XtX   * sl) - as.numeric(bl %*% XtXbl)
    # trXtXsqM <- trXtXsqM + sum(diag_XtXsq * sl) - as.numeric(bl %*% XtXsqbl)
    trXtXM   <- trXtXM   + sum(diag_XtX   * sl) - sum(bl * as.numeric(XtXbl))
    trXtXsqM <- trXtXsqM + sum(diag_XtXsq * sl) - sum(bl * as.numeric(XtXsqbl))
  }
  
  # x0 = yty - 2 b'Xty + tr(X'X M)
  # x1 = sum(Xty^2) - 2 (X'Xb)'Xty + tr(X'X^2 M)
  # x0 <- as.numeric(yty) - 2 * as.numeric(b %*% Xty) + trXtXM
  # x1 <- sum(Xty^2)      - 2 * as.numeric(XtXb %*% Xty) + trXtXsqM
  # x1 <- sum(Xty^2) - 2 * as.numeric(crossprod(XtXb, Xty)) + trXtXsqM
  x0 <- as.numeric(yty) - 2 * sum(b * as.numeric(Xty)) + trXtXM
  x1 <- sum(Xty^2)      - 2 * sum(as.numeric(XtXb) * as.numeric(Xty)) + trXtXsqM
  
  x <- c(x0, x1)
  
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]; tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n; tausq <- 0
    }
    if (verbose) cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

MoMK <- function(PIP, mu, omega, sigmasq, tausq, n,
                 Z_new, D_Z, y_new,          # replaces XtX, XtXsq
                 Xty, yty,
                 diagXtX, diagXtXsq,         # precomputed diagonals
                 trK, trXtX, trXtKX,
                 est_sigmasq, est_tausq, verbose) {
  
  L <- ncol(mu)
  
  A <- matrix(c(n,     trK,
                trXtX, trXtKX), 2, 2, byrow = TRUE)
  
  b  <- rowSums(mu * PIP)
  Zb <- as.numeric(Z_new %*% b)           # n-vector, O(np)
  
  # b'(X'X)b = ||Zb||^2
  # b'(X'X)^2 b = (Zb)' D_Z (Zb) = sum(D_Z * Zb^2)
  trXtXM   <- sum(Zb^2)
  trXtXsqM <- sum(D_Z * Zb^2)
  
  for (l in seq_len(L)) {
    bl  <- mu[, l] * PIP[, l]
    sl  <- PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
    Zbl <- as.numeric(Z_new %*% bl)       # n-vector, O(np)
    
    trXtXM   <- trXtXM   + sum(diagXtX   * sl) - sum(Zbl^2)
    trXtXsqM <- trXtXsqM + sum(diagXtXsq * sl) - sum(D_Z * Zbl^2)
  }
  
  x0 <- yty       - 2 * sum(b * Xty)              + trXtXM
  # sum(XtXb * Xty) = b'(X'X)(X'y) = (Zb)' D_Z y_new
  x1 <- sum(Xty^2) - 2 * sum(D_Z * Zb * y_new)   + trXtXsqM
  
  x <- c(x0, x1)
  
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]; tausq <- sol[2]
    } else {
      sigmasq <- x[1] / n; tausq <- 0
    }
    if (verbose) cat(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigmasq, tausq))
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

# ──────────────────────────────────────────────────────────────────────────────
# MoMKv2
#
# Method-of-moments estimator for (sigma^2, tau^2) using weight matrices
#   W1 = I   =>  E[(y-Xb)'I(y-Xb)]   = sigma^2 * n       + tau^2 * tr(K)
#   W2 = K   =>  E[(y-Xb)'K(y-Xb)]   = sigma^2 * tr(K)   + tau^2 * tr(K^2)
#
# All computations are done in the eigenbasis (Z = Q'X, y_tilde = Q'y),
# which makes the K-weighted quantities cheap: multiplying by K in the
# original basis is just multiplying by diag(eigvals) in the eigenbasis.
#
# Arguments
#   PIP, mu, omega  : p x L variational parameters
#   sigmasq, tausq  : current variance estimates (returned unchanged if not updated)
#   n               : sample size
#   Z               : n x p matrix Q'X  (eigenbasis-projected genotypes)
#   eigvals         : length-n eigenvalues of K (same ordering as rows of Z)
#   y_tilde         : length-n vector Q'y
#   Xty             : length-p vector X'y  (= Z'y_tilde, precomputed in susieKv2)
#   XtKy            : length-p vector X'Ky (= Z'(eigvals * y_tilde))
#   yty             : scalar y'y
#   ytKy            : scalar y'Ky (= sum(eigvals * y_tilde^2))
#   diagXtX         : length-p diagonal of X'X  (= colSums(Z^2))
#   diagXtKX        : length-p diagonal of X'KX (= Z^2' eigvals, i.e. crossprod(Z^2, eigvals))
#   trK             : tr(K) = sum(eigvals)
#   trK2            : tr(K^2) = sum(eigvals^2)
#   est_sigmasq     : logical
#   est_tausq       : logical
#   verbose         : logical
# ──────────────────────────────────────────────────────────────────────────────
MoMKv2 <- function(PIP, mu, omega, sigmasq, tausq, n,
                   Z, eigvals, y_tilde,
                   Xty, XtKy, yty, ytKy,
                   diagXtX, diagXtKX,
                   trK, trK2,
                   est_sigmasq, est_tausq, verbose) {
  
  L <- ncol(mu)
  
  # ── moment matrix A ──────────────────────────────────────────────────────────
  # A[i,j] = tr(Wi Cj) where C1 = I, C2 = K
  #   row 1 (W1=I):  [tr(I),    tr(K)  ] = [n,   trK ]
  #   row 2 (W2=K):  [tr(K),    tr(K^2)] = [trK, trK2]
  A <- matrix(c(n,   trK,
                trK, trK2), 2, 2, byrow = TRUE)
  
  # ── posterior mean vector b and its projections ──────────────────────────────
  b   <- rowSums(mu * PIP)            # length-p
  Zb  <- as.numeric(Z %*% b)         # length-n  (X*beta in eigenbasis)
  
  # ── tr(X'X M) and tr(X'KX M) ─────────────────────────────────────────────────
  # Second moment matrix M = b b' + sum_l diag(s_l - b_l o b_l)
  # tr(X'X   M) = b'X'X b   + sum_l [ diag(X'X)   . s_l - ||Z b_l||^2        ]
  # tr(X'KX  M) = b'X'KX b  + sum_l [ diag(X'KX)  . s_l - sum_j lj (Zbl_j)^2 ]
  trXtXM  <- sum(Zb^2)                   # b'X'Xb  = ||Zb||^2
  trXtKXM <- sum(eigvals * Zb^2)         # b'X'KXb = sum_j lj (Zb_j)^2
  
  for (l in seq_len(L)) {
    bl  <- mu[, l] * PIP[, l]                             # p-vector
    sl  <- PIP[, l] * (mu[, l]^2 + 1 / omega[, l])       # p-vector (2nd moment)
    Zbl <- as.numeric(Z %*% bl)                           # n-vector
    
    trXtXM  <- trXtXM  + sum(diagXtX  * sl) - sum(Zbl^2)
    trXtKXM <- trXtKXM + sum(diagXtKX * sl) - sum(eigvals * Zbl^2)
  }
  
  # ── residual quadratic forms (observed moments) ───────────────────────────────
  # x0 = E[(y-Xb)'I(y-Xb)]  = y'y - 2 b'X'y + tr(X'X M)
  # x1 = E[(y-Xb)'K(y-Xb)]  = y'Ky - 2 b'X'Ky + tr(X'KX M)
  x0 <- yty  - 2 * sum(b * Xty)  + trXtXM
  x1 <- ytKy - 2 * sum(b * XtKy) + trXtKXM
  
  x <- c(x0, x1)
  
  # ── solve ─────────────────────────────────────────────────────────────────────
  if (est_tausq) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigmasq <- sol[1]
      tausq   <- sol[2]
    } else {
      # fallback: pure noise, no structured component
      sigmasq <- x[1] / n
      tausq   <- 0
    }
    if (verbose) cat(sprintf("  Update (sigma^2, tau^2) to (%f, %e)\n", sigmasq, tausq))
  } else if (est_sigmasq) {
    sigmasq <- (x[1] - trK * tausq) / n
    if (verbose) cat(sprintf("  Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

# ──────────────────────────────────────────────────────────────────────────────
# VCupdateK
#
# Estimates (sigma^2, tau^2) by maximising the ELBO contribution from the
# likelihood with respect to the variance components, given the current
# variational posterior q.
#
# In the eigenbasis of K the ELBO decomposes as:
#
#   ell(sigma^2, tau^2) = -1/2 * sum_i [ log(c_i) + D_i / c_i ]
#
# where c_i = tau^2 * lambda_i + sigma^2  and
#       D_i = E_q[(y_tilde_i - (Z beta)_i)^2]
#           = r_i^2 + (Z Sigma_q Z')_ii
#           = (eigenbasis residual)^2 + (eigenbasis posterior variance)_i
#
# The stationarity conditions are:
#   sum_i w_i^2 * D_i = sum_i w_i           (d/d sigma^2 = 0)
#   sum_i lam_i * w_i^2 * D_i = sum_i lam_i * w_i  (d/d tau^2 = 0)
#
# where w_i = 1/c_i. These are solved by IRLS: at each step, holding w_i
# fixed at current values, the conditions become a 2x2 linear system whose
# solution is the next (sigma^2, tau^2). The weights are then updated and
# the process repeats.
#
# This is strictly better than MoM (W=I, W=K) for general K because it
# weights each eigendirection by 1/c_i^2 (Fisher information), suppressing
# the large-eigenvalue directions that dominate unweighted MoM but carry
# little information about the variance components.
#
# Arguments
#   PIP, mu, omega   : p x L variational parameters
#   sigmasq, tausq   : current variance estimates (used as starting point)
#   Z                : n x p matrix Q'X
#   eigvals          : length-n eigenvalues of K (ascending order, matching Z)
#   y_tilde          : length-n vector Q'y
#   est_sigmasq      : logical — whether to update sigma^2
#   est_tausq        : logical — whether to update tau^2
#   maxiter_vc, tol_vc : IRLS convergence controls
# ──────────────────────────────────────────────────────────────────────────────
VCupdateK <- function(PIP, mu, omega, sigmasq, tausq,
                      Z, eigvals, y_tilde,
                      est_sigmasq, est_tausq,
                      maxiter_vc = 20, tol_vc = 1e-8, verbose) {
  
  L <- ncol(mu)
  n <- length(y_tilde)
  
  # ── compute D_i = E_q[(y_tilde_i - (Z beta)_i)^2] ──────────────────────────
  # D = r^2 + diag(Z Sigma_q Z')
  # r_i = y_tilde_i - (Zb)_i    (eigenbasis residual)
  # diag(Z Sigma_q Z')_i = sum_l [ (Z^2 s_l)_i - (Z b_l)_i^2 ]
  b <- rowSums(mu * PIP)
  r <- y_tilde - as.numeric(Z %*% b)
  
  v <- numeric(n)
  for (l in seq_len(L)) {
    bl  <- mu[, l] * PIP[, l]
    sl  <- PIP[, l] * (mu[, l]^2 + 1 / omega[, l])
    Zbl <- as.numeric(Z %*% bl)
    v   <- v + as.numeric(Z^2 %*% sl) - Zbl^2
  }
  D <- r^2 + v
  
  # ── IRLS ────────────────────────────────────────────────────────────────────
  for (iter in seq_len(maxiter_vc)) {
    c_i  <- tausq * eigvals + sigmasq  # model variance per eigendirection
    w2   <- 1 / c_i^2                 # Fisher information weights
    lw2  <- eigvals * w2
    l2w2 <- eigvals * lw2
    
    # 2x2 WLS system:  B [sigma^2, tau^2]' = d
    # B = X'WX where W = diag(w2), design matrix = [1, lambda_i]
    B <- matrix(c(sum(w2),   sum(lw2),
                  sum(lw2),  sum(l2w2)), 2, 2, byrow = TRUE)
    d <- c(sum(w2 * D), sum(lw2 * D))
    
    if (!est_tausq) {
      sig_new <- (d[1] - B[1, 2] * tausq) / B[1, 1]
      tau_new <- tausq
    } else {
      sol <- tryCatch(solve(B, d), error = function(e) c(NA_real_, NA_real_))
      if (anyNA(sol) || sol[1] <= 0 || sol[2] < 0) {
        # boundary: tau^2 = 0, sigma^2 = weighted mean of D
        sig_new <- d[1] / B[1, 1]
        tau_new <- 0
      } else {
        sig_new <- sol[1]
        tau_new <- sol[2]
      }
    }
    
    delta   <- max(abs(sig_new - sigmasq), abs(tau_new - tausq))
    sigmasq <- sig_new
    tausq   <- tau_new
    
    if (verbose) cat(sprintf("  VC iter %d: sigma^2=%f tau^2=%e (delta=%e)\n",
                             iter, sigmasq, tausq, delta))
    if (delta < tol_vc) break
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}
