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
MoMK <- function(PIP, mu, omega, sigmasq, tausq, n, eigvals,
                 XtX_diag, XTKX_diag, Xty, XTKy, yty, yKy, trK, trK2,
                 est_sigmasq, est_tausq, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  # Use byrow=TRUE to match Python row order
  A <- matrix(c(n, trK, trK, trK2), 2, 2, byrow = TRUE)
  
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
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
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

MoMK <- function(PIP, mu, omega, sigmasq, tausq, n, eigvals,
                 XtX_diag, XTKX_diag, Xty, XTKy, yty, yKy, trK, trK2,
                 est_sigmasq, est_tausq, verbose) {
  
  p <- nrow(mu); L <- ncol(mu)
  
  # Use byrow=TRUE to match Python row order
  A <- matrix(c(n, trK, trK, trK2), 2, 2, byrow = TRUE)
  
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
    sigmasq <- (x[1] - A[1, 2] * tausq) / n
    if (verbose) cat(sprintf("Update sigma^2 to %f\n", sigmasq))
  }
  
  list(sigmasq = sigmasq, tausq = tausq)
}

