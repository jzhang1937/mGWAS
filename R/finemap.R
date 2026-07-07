# run susie on each
#' Run SuSiE on each locus using individual-level data
#'
#' @param loci     output of define_loci()
#' @param X        n x p genotype matrix
#' @param y        length-n phenotype vector
#' @param L        max number of causal signals per locus (default 10)
#' @param min_purity  min purity to retain a credible set (default 0.5)
#'
#' @return list of results, one per locus, each with:
#'   $fit  : raw susie fit object
#'   $cs   : credible sets with original variant indices
#'   $pip  : named vector of PIPs for locus variants (original indices)

run_susie_loci <- function(loci, X, y, L = 10, min_purity = 0.5) {
  
  if (!requireNamespace("susieR", quietly = TRUE))
    stop("Please install susieR: install.packages('susieR')")
  
  results <- vector("list", length(loci))
  
  for (k in seq_along(loci)) {
    
    idx <- loci[[k]]$variants
    message(sprintf("\nLocus %d: %d variants, lead(s): %s",
                    k, length(idx), paste(loci[[k]]$leads, collapse = ", ")))
    
    fit <- tryCatch(
      susieR::susie(
        X = X[, idx, drop = FALSE],
        y = y,
        L = L
      ),
      error = function(e) {
        message(sprintf("  SuSiE failed: %s", e$message))
        NULL
      }
    )
    
    if (is.null(fit)) {
      results[[k]] <- list(fit = NULL, cs = NULL, pip = NULL)
      next
    }
    
    # Remap credible set indices from locus-local to original
    cs_original <- NULL
    if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
      
      purity    <- fit$sets$purity$min.abs.corr
      keep      <- which(purity >= min_purity)
      n_dropped <- length(fit$sets$cs) - length(keep)
      
      if (n_dropped > 0)
        message(sprintf("  Dropped %d low-purity CS.", n_dropped))
      
      cs_original <- lapply(fit$sets$cs[keep], function(cs_idx) idx[cs_idx])
      message(sprintf("  %d credible set(s) retained.", length(cs_original)))
      
    } else {
      message("  No credible sets found.")
    }
    
    pip_original <- setNames(fit$pip, idx)
    
    results[[k]] <- list(
      fit = fit,
      cs  = cs_original,
      pip = pip_original
    )
  }
  
  results
}

#' Coordinate-free locus definition by greedy clumping + neighborhood expansion
#'
#' @param pvals    length-p numeric vector of p-values
#' @param R2       p x p matrix of r^2 values between variants
#' @param sig_thr  significance threshold for lead variants (default 5e-8)
#' @param r2_thr   r^2 threshold for clumping and neighborhood (default 0.1)
#'
#' @return list of loci, each with:
#'   $lead     : index of lead variant
#'   $variants : all variant indices in locus (lead + correlated neighbors)

define_loci <- function(pvals, R2, sig_thr = 5e-8, r2_thr = 0.1) {
  
  p <- length(pvals)
  stopifnot(nrow(R2) == p, ncol(R2) == p)
  
  # --- 1. Greedy clumping to find lead variants ---
  # Sort significant variants by p-value (best first)
  candidates  <- order(pvals)
  candidates  <- candidates[pvals[candidates] < sig_thr]
  
  if (length(candidates) == 0) stop("No variants pass significance threshold.")
  message(sprintf("%d variants pass significance threshold.", length(candidates)))
  
  leads      <- c()
  clumped_out <- c()
  
  for (i in candidates) {
    if (i %in% clumped_out) next
    leads       <- c(leads, i)
    correlated  <- which(R2[i, ] > r2_thr)
    clumped_out <- union(clumped_out, correlated)
  }
  
  message(sprintf("%d lead variant(s) after clumping.", length(leads)))
  
  # --- 2. Expand each lead to its neighborhood ---
  # neighborhood: all variants with r^2 > r2_thr to the lead
  neighborhoods <- lapply(leads, function(lead) which(R2[lead, ] > r2_thr))
  
  # --- 3. Merge neighborhoods whose leads are correlated ---
  # Build membership vector over leads; merge if any two leads share r^2 > r2_thr
  n_leads    <- length(leads)
  membership <- seq_len(n_leads)
  
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (i in seq_len(n_leads)) {
      for (j in seq_len(n_leads)) {
        if (i >= j) next
        if (R2[leads[i], leads[j]] > r2_thr && membership[i] != membership[j]) {
          membership[membership == membership[j]] <- membership[i]
          changed <- TRUE
        }
      }
    }
  }
  
  # Relabel components 1..K
  comp_ids   <- unique(membership)
  n_loci     <- length(comp_ids)
  message(sprintf("%d locus/loci after merging correlated leads.", n_loci))
  
  # --- 4. Build loci: union of neighborhoods within each component ---
  loci <- vector("list", n_loci)
  
  for (k in seq_along(comp_ids)) {
    in_comp      <- which(membership == comp_ids[k])
    comp_leads   <- leads[in_comp]
    comp_vars    <- sort(unique(unlist(neighborhoods[in_comp])))
    
    loci[[k]] <- list(
      leads    = comp_leads,
      variants = comp_vars
    )
    
    message(sprintf("  Locus %d: %d lead(s), %d variants",
                    k, length(comp_leads), length(comp_vars)))
  }
  
  loci
}

#' Coordinate-free locus definition by greedy clumping + neighborhood expansion
#' (memory-light version: computes r^2 on the fly from genotypes, never
#' materializes a p x p matrix)
#'
#' @param pvals    length-p numeric vector of p-values
#' @param X        n x p genotype matrix (samples x variants), matching
#'                 column order of pvals. Must be complete (no NA).
#' @param sig_thr  significance threshold for lead variants (default 5e-8)
#' @param r2_thr   r^2 threshold for clumping and neighborhood (default 0.1)
#' @param already_scaled  set TRUE if X is already column-standardized
#'                 (mean 0, sd 1) to skip re-scaling it internally
#'
#' @return list of loci, each with:
#'   $lead     : index of lead variant
#'   $variants : all variant indices in locus (lead + correlated neighbors)
define_loci_fast <- function(pvals, X, sig_thr = 5e-8, r2_thr = 0.1,
                             already_scaled = FALSE) {
  
  p <- length(pvals)
  n <- nrow(X)
  stopifnot(ncol(X) == p)
  
  # --- check names/colnames align ---
  if (is.null(names(pvals)) || is.null(colnames(X))) {
    stop("Both pvals (names) and X (colnames) must be named to check alignment.")
  }
  if (!identical(names(pvals), colnames(X))) {
    stop("names(pvals) and colnames(X) do not match (order or content differs). ",
         "Reorder/rename before calling define_loci_fast().")
  }
  
  # --- 0. Standardize once: O(n*p), same memory order as X itself ---
  # After scaling, cor(Xs[,i], Xs[,j]) = crossprod(Xs[,i], Xs[,j]) / (n - 1)
  Xs <- if (already_scaled) X else scale(X)
  
  # --- 1. Greedy clumping to find lead variants ---
  candidates <- order(pvals)
  candidates <- candidates[pvals[candidates] < sig_thr]
  
  if (length(candidates) == 0) stop("No variants pass significance threshold.")
  message(sprintf("%d variants pass significance threshold.", length(candidates)))
  
  leads       <- integer(0)
  clumped_out <- logical(p)   # O(p) bit vector instead of union() on growing vectors
  
  for (i in candidates) {
    if (clumped_out[i]) next
    leads <- c(leads, i)
    clumped_out[i] <- TRUE
    
    # r^2 of this lead against remaining CANDIDATES only (not all p) 
    still_open <- candidates[!clumped_out[candidates]]
    if (length(still_open) > 0) {
      r <- as.vector(crossprod(Xs[, i], Xs[, still_open, drop = FALSE])) / (n - 1)
      correlated <- still_open[r^2 > r2_thr]
      clumped_out[correlated] <- TRUE
    }
  }
  
  message(sprintf("%d lead variant(s) after clumping.", length(leads)))
  
  # --- 2. Expand each lead to its neighborhood across ALL p variants ---
  neighborhoods <- lapply(leads, function(lead) {
    r <- as.vector(crossprod(Xs[, lead], Xs)) / (n - 1)
    which(r^2 > r2_thr)
  })
  
  # --- 3. Merge neighborhoods whose leads are correlated ---
  K <- length(leads)
  if (K > 1) {
    R2_leads <- (crossprod(Xs[, leads, drop = FALSE]) / (n - 1))^2
  } else {
    R2_leads <- matrix(1, 1, 1)
  }
  
  membership <- seq_len(K)
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (i in seq_len(K)) {
      for (j in seq_len(K)) {
        if (i >= j) next
        if (R2_leads[i, j] > r2_thr && membership[i] != membership[j]) {
          membership[membership == membership[j]] <- membership[i]
          changed <- TRUE
        }
      }
    }
  }
  
  comp_ids <- unique(membership)
  n_loci   <- length(comp_ids)
  message(sprintf("%d locus/loci after merging correlated leads.", n_loci))
  
  # --- 4. Build loci: union of neighborhoods within each component ---
  loci <- vector("list", n_loci)
  
  for (k in seq_along(comp_ids)) {
    in_comp    <- which(membership == comp_ids[k])
    comp_leads <- leads[in_comp]
    comp_vars  <- sort(unique(unlist(neighborhoods[in_comp])))
    
    loci[[k]] <- list(
      leads    = comp_leads,
      variants = comp_vars
    )
    
    message(sprintf("  Locus %d: %d lead(s), %d variants",
                    k, length(comp_leads), length(comp_vars)))
  }
  
  loci
}

#' Combine SuSiE results across loci into a single tidy data frame
#'
#' @param results     output of run_susie_loci()
#' @param loci        output of define_loci()
#' @param var_names   length-p character vector of variant names
#'
#' @return data frame with one row per variant that appears in any locus, with columns:
#'   variant_idx, variant_name, locus, is_lead, pip, cs (NA if not in a CS)

compile_results <- function(results, loci, var_names) {
  
  rows <- lapply(seq_along(loci), function(k) {
    
    idx     <- loci[[k]]$variants
    leads   <- loci[[k]]$leads
    pip     <- results[[k]]$pip        # named by original index
    cs_list <- results[[k]]$cs         # list of integer vectors (original indices)
    
    # Map each variant to its CS (NA if not in any CS)
    cs_membership <- rep(NA_integer_, length(idx))
    if (!is.null(cs_list) && length(cs_list) > 0) {
      for (cs_k in seq_along(cs_list)) {
        in_cs <- idx %in% cs_list[[cs_k]]
        cs_membership[in_cs] <- cs_k
      }
    }
    
    data.frame(
      variant_idx  = idx,
      variant_name = var_names[idx],
      locus        = k,
      is_lead      = idx %in% leads,
      pip          = if (is.null(pip)) NA_real_ else as.numeric(pip),
      cs           = cs_membership,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, rows)
}

#' Run susieKv2 on each locus
#'
#' @param loci       output of define_loci()
#' @param data       list with $X (n x p), $y (length-n), $tree (phylo object or NULL)
#' @param L          max number of causal signals per locus (default 10)
#' @param min_pip    PIP threshold to define a credible set (variants with PIP > threshold)
#' @param ...        additional arguments passed to susieKv2
#'
#' @return list of results, one per locus, each with:
#'   $fit         : raw susieKv2 fit object
#'   $cs          : credible sets with original variant indices
#'   $pip         : named vector of marginal PIPs for locus variants (original indices)

run_susieKv2_loci <- function(loci, data, L = 10, min_pip = 0.95, ...) {
  
  results <- vector("list", length(loci))
  
  for (k in seq_along(loci)) {
    
    idx <- loci[[k]]$variants
    message(sprintf("\nLocus %d: %d variants, lead(s): %s",
                    k, length(idx), paste(loci[[k]]$leads, collapse = ", ")))
    
    # Subset data to locus variants
    locus_data <- list(
      X    = data$X[, idx, drop = FALSE],
      y    = data$y,
      tree = data$tree
    )
    
    fit <- tryCatch(
      mGWAS::susieKv2(data = locus_data, L = L, ...),
      error = function(e) {
        message(sprintf("  susieKv2 failed: %s", e$message))
        NULL
      }
    )
    
    if (is.null(fit)) {
      results[[k]] <- list(fit = NULL, cs = NULL, pip = NULL)
      next
    }
    
    # susieKv2 returns marginalPIP (length-p vector, named by colnames of locus X)
    pip <- fit$marginalPIP
    
    # Credible sets: for each of the L components, collect variants above min_pip
    cs_original <- vector("list", L)
    for (l in seq_len(L)) {
      in_cs <- which(fit$PIP[, l] > min_pip)
      if (length(in_cs) > 0) {
        cs_original[[l]] <- idx[in_cs]
      }
    }
    # Drop empty components
    cs_original <- Filter(Negate(is.null), cs_original)
    
    if (length(cs_original) > 0) {
      message(sprintf("  %d credible set(s) found.", length(cs_original)))
    } else {
      message("  No credible sets found.")
      cs_original <- NULL
    }
    
    results[[k]] <- list(
      fit = fit,
      cs  = cs_original,
      pip = setNames(pip, idx)
    )
  }
  
  results
}

#' Define loci by hierarchical clustering on r² distance (no positional structure)
#'
#' @param R2       p x p matrix of r² values between variants
#' @param r2_thr   r² threshold — variants with r² > r2_thr end up in same cluster (default 0.3)
#' @param method   linkage method: "average", "complete", or "single" (default "average")
#'
#' @return list of loci, each with:
#'   $leads    : index of the variant with highest mean r² to cluster members (medoid)
#'   $variants : all variant indices in locus

define_loci_clust <- function(R2, r2_thr = 0.3, method = "average") {
  
  p <- nrow(R2)
  
  D  <- 1 - R2
  D[D < 0] <- 0
  diag(D)  <- 0
  
  hc     <- hclust(as.dist(D), method = method)
  labels <- cutree(hc, h = 1 - r2_thr)
  
  n_loci <- max(labels)
  message(sprintf("%d loci defined from %d variants (r2_thr = %.2f, linkage = %s)",
                  n_loci, p, r2_thr, method))
  
  loci <- vector("list", n_loci)
  
  for (k in seq_len(n_loci)) {
    
    idx <- which(labels == k)
    
    # medoid: variant with highest mean r² to all other cluster members
    if (length(idx) == 1) {
      lead <- idx
    } else {
      mean_r2 <- rowMeans(R2[idx, idx, drop = FALSE])
      lead    <- idx[which.max(mean_r2)]
    }
    
    loci[[k]] <- list(
      leads    = lead,
      variants = idx
    )
  }
  
  cluster_sizes <- sapply(loci, function(l) length(l$variants))
  
  message(sprintf("  Cluster sizes: min=%d, median=%d, max=%d",
                  min(cluster_sizes),
                  as.integer(median(cluster_sizes)),
                  max(cluster_sizes)))
  
  loci
}

#' Wrapper to run locus definition + SuSiE and compile into a standard output
#'
#' @param X              n x p genotype matrix
#' @param y              length-n phenotype vector
#' @param define_fn      locus definition function, must return list with $variants and $leads
#' @param susie_fn       susie runner function, must accept (loci, X, y, ...)
#' @param var_names      optional length-p character vector of variant names
#' @param ...            additional arguments passed to susie_fn
#'
#' @return list with:
#'   $pip      : length-p vector of PIPs (max across loci if overlapping)
#'   $pip_mean : length-p vector of mean PIPs across loci
#'   $pip_all  : length-p list of all PIP estimates per variant
#'   $pip_n    : length-p integer vector of how many loci each variant appeared in
#'   $cs       : list of credible sets as original variant indices
#'   $loci     : raw loci definition
#'   $results  : raw per-locus susie results

run_susie_clustered <- function(X, y, define_fn, susie_fn, var_names = NULL, ...) {
  
  p <- ncol(X)
  if (is.null(var_names)) var_names <- as.character(seq_len(p))
  
  # --- 1. Define loci ---
  message("=== Defining loci ===")
  loci <- define_fn()
  
  # --- 2. Run SuSiE on each locus ---
  message("\n=== Running SuSiE per locus ===")
  results <- susie_fn(loci, X, y, ...)
  
  # --- 3. Compile PIPs, handling overlapping loci ---
  message("\n=== Compiling results ===")
  
  pip_collect <- vector("list", p)
  names(pip_collect) <- var_names
  
  for (k in seq_along(loci)) {
    idx   <- loci[[k]]$variants
    pip_k <- results[[k]]$pip
    if (!is.null(pip_k)) {
      for (i in seq_along(idx)) {
        pip_collect[[idx[i]]] <- c(pip_collect[[idx[i]]], as.numeric(pip_k)[i])
      }
    }
  }
  
  pip_max  <- sapply(pip_collect, function(x) if (is.null(x)) NA_real_ else max(x))
  pip_mean <- sapply(pip_collect, function(x) if (is.null(x)) NA_real_ else mean(x))
  pip_n    <- sapply(pip_collect, function(x) if (is.null(x)) 0L else length(x))
  
  n_repeated <- sum(pip_n > 1)
  if (n_repeated > 0)
    message(sprintf("  %d variant(s) appeared in multiple loci.", n_repeated))
  
  # --- 4. Compile credible sets globally ---
  cs_global  <- list()
  cs_counter <- 1
  
  for (k in seq_along(loci)) {
    cs_k <- results[[k]]$cs
    if (!is.null(cs_k) && length(cs_k) > 0) {
      for (cs in cs_k) {
        cs_global[[cs_counter]] <- cs
        cs_counter <- cs_counter + 1
      }
    }
  }
  
  n_cs <- length(cs_global)
  message(sprintf("Total credible sets across all loci: %d", n_cs))
  message(sprintf("Variants with PIP > 0.5: %d", sum(pip_max > 0.5, na.rm = TRUE)))
  message(sprintf("Variants with PIP > 0.9: %d", sum(pip_max > 0.9, na.rm = TRUE)))
  
  list(
    pip      = pip_max,
    pip_mean = pip_mean,
    pip_all  = pip_collect,
    pip_n    = pip_n,
    cs       = cs_global,
    loci     = loci,
    results  = results
  )
}

#' Run SuSiE on cluster representatives (medoids) only
#'
#' Reduces dimensionality by clustering correlated variants, picking one
#' representative per cluster, and running a single SuSiE fit on just the
#' representatives. Useful as a coarse screening step before full fine-mapping.
#'
#' @param X          n x p genotype matrix
#' @param y          length-n phenotype vector
#' @param define_fn  clustering function (zero-arg closure), must return list
#'                    of loci each with $variants and $leads (as in define_loci_clust)
#' @param L          max number of causal signals (default 10)
#' @param min_purity min purity to retain a credible set (default 0.5)
#' @param var_names  optional length-p character vector of variant names
#' @param ...         additional arguments passed to susieR::susie
#'
#' @return list with:
#'   $pip        : length-p vector, NA for non-representatives, PIP for representatives
#'   $pip_rep    : named vector of PIPs, only for representative variants
#'   $cs         : credible sets, as original variant indices
#'   $fit        : raw susie fit object (on representatives only)
#'   $loci       : raw cluster definitions
#'   $rep_idx    : indices of variants used as representatives

run_susie_representatives <- function(X, y, define_fn, L = 10,
                                      min_purity = 0.5, var_names = NULL, ...) {
  
  if (!requireNamespace("susieR", quietly = TRUE))
    stop("Please install susieR: install.packages('susieR')")
  
  p <- ncol(X)
  if (is.null(var_names)) var_names <- as.character(seq_len(p))
  
  # --- 1. Define clusters ---
  message("=== Defining clusters ===")
  loci <- define_fn()
  n_clusters <- length(loci)
  
  # --- 2. Extract one representative (lead/medoid) per cluster ---
  rep_idx <- sapply(loci, function(l) l$leads[1])  # take first if multiple leads
  
  if (any(duplicated(rep_idx)))
    warning("Duplicate representatives detected across clusters.")
  
  message(sprintf("%d clusters -> %d representative variants (dimensionality %d -> %d)",
                  n_clusters, length(rep_idx), p, length(rep_idx)))
  
  # --- 3. Run a single SuSiE fit on just the representatives ---
  message("\n=== Running SuSiE on representatives ===")
  
  fit <- tryCatch(
    susieR::susie(
      X = X[, rep_idx, drop = FALSE],
      y = y,
      L = L,
      ...
    ),
    error = function(e) {
      message(sprintf("  SuSiE failed: %s", e$message))
      NULL
    }
  )
  
  if (is.null(fit)) {
    return(list(pip = rep(NA_real_, p), pip_rep = NULL, cs = NULL,
                fit = NULL, loci = loci, rep_idx = rep_idx))
  }
  
  # --- 4. Remap credible sets from representative-local to original indices ---
  cs_original <- NULL
  if (!is.null(fit$sets$cs) && length(fit$sets$cs) > 0) {
    
    purity    <- fit$sets$purity$min.abs.corr
    keep      <- which(purity >= min_purity)
    n_dropped <- length(fit$sets$cs) - length(keep)
    
    if (n_dropped > 0)
      message(sprintf("  Dropped %d low-purity CS.", n_dropped))
    
    cs_original <- lapply(fit$sets$cs[keep], function(cs_idx) rep_idx[cs_idx])
    message(sprintf("  %d credible set(s) retained.", length(cs_original)))
    
  } else {
    message("  No credible sets found.")
  }
  
  # --- 5. Build length-p PIP vector (NA for non-representatives) ---
  pip_rep <- setNames(fit$pip, rep_idx)
  
  pip_global <- rep(NA_real_, p)
  names(pip_global) <- var_names
  pip_global[rep_idx] <- fit$pip
  
  message(sprintf("Representatives with PIP > 0.5: %d", sum(fit$pip > 0.5)))
  message(sprintf("Representatives with PIP > 0.9: %d", sum(fit$pip > 0.9)))
  
  list(
    pip     = pip_global,
    pip_rep = pip_rep,
    cs      = cs_original,
    fit     = fit,
    loci    = loci,
    rep_idx = rep_idx
  )
}