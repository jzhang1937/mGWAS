fdp <- function(selected, nonnulls) {
  false_positives <- setdiff(selected, nonnulls)
  n_selected <- length(selected)
  if (n_selected == 0) {
    return(0)
  } else {
    return(length(false_positives) / n_selected)
  }
}

power <- function(selected, nonnulls) {
  false_negatives <- setdiff(nonnulls, selected)
  n_nonnulls <- length(nonnulls)
  return((n_nonnulls - length(false_negatives)) / n_nonnulls)
}

#' Title
#'
#' @param tree a phylogenetic tree
#' @param labels labels of the tree to keep
#'
#' @returns a pruned tree which only keeps data points in labels
#' @export
#'
#' @examples
tree_prune <- function(tree, labels) {
  tree_pruned <- ape::keep.tip(tree, labels)
  tree_pruned
}

#' Title
#'
#' @param treeWAS_object a treeWAS output object
#'
#' @returns list object with each scores p-values
#' @export
#'
#' @examples
process_treeWAS <- function(treeWAS_object) {
  scores <- c("terminal", "simultaneous", "subsequent")
  results <- vector(mode = "list", length(scores))
  names(results) <- scores
  for (score in scores) {
    tree_result <- treeWAS_object[[score]]
    # Sort absolute values
    sorted <- sort(abs(tree_result$corr.sim))
    eps <- .Machine$double.eps 
    # find tail
    counts <- length(sorted) - findInterval(abs(tree_result$corr.dat) - eps, sorted)
    p.vals <- (counts + 1) / (length(sorted) + 1)
    names(p.vals) <- names(tree_result$corr.dat)
    results[[score]] <- p.vals
  }
  results
}

#' Title
#'
#' @param treeWAS_object a treeWAS output object
#'
#' @returns list object with each scores p-values directly from treeWAS
#' @export
#'
#' @examples
raw_treeWAS <- function(treeWAS_object) {
  scores <- c("terminal", "simultaneous", "subsequent")
  results <- vector(mode = "list", length(scores))
  names(results) <- scores
  for (score in scores) {
    tree_result <- treeWAS_object[[score]]
    p.vals <- tree_result$p.vals
    results[[score]] <- p.vals
  }
  results
}



#' Title
#'
#' @param X 
#' @param S.hat 
#' @param S.star 
#'
#' @returns
#' @export
#'
#' @examples
compute_subspace_tp <- function(X, S.hat, S.star) {
  if (length(S.hat) == 0 | length(S.star) == 0) {
    return(0)
  } 
  X.S.hat <- X[,S.hat]
  X.S.star <- X[,S.star]
  A <- X.S.hat %*% MASS::ginv(t(X.S.hat) %*% X.S.hat) %*% t(X.S.hat) %*% 
    X.S.star %*% MASS::ginv(t(X.S.star) %*% X.S.star) %*% t(X.S.star)
  sum(diag(A))
}



#' Flag variants whose marginal PIP changed
#'
#' @param PIP matrix of PIPs
#' @param PIP_prev matrix of prev PIPs
#' @param threshold minimum absolute change in marginal PIP to flag
#' @param variant_names optional vector of variant names, same length/order as nrow(PIP)
#'
#' @returns list with changed_idx, changed_variants, and marginal_diff
#' @export
#'
#' @examples
check_pip_changes_marginal <- function(PIP, PIP_prev, threshold = 1e-4, variant_names = NULL) {
  marginal_diff <- (1 - apply(1 - PIP, 1, prod)) - (1 - apply(1 - PIP_prev, 1, prod))
  idx <- which(abs(marginal_diff) > threshold)
  
  list(
    changed_idx      = idx,
    changed_variants = if (!is.null(variant_names)) variant_names[idx] else idx,
    marginal_diff    = marginal_diff
  )
}

#' Flag variants whose PIP changed in any of the L effects
#'
#' @param PIP matrix of PIPs
#' @param PIP_prev matrix of prev PIPs
#' @param threshold minimum absolute change, in any single effect, to flag
#' @param variant_names optional vector of variant names, same length/order as nrow(PIP)
#'
#' @returns list with changed_idx, changed_variants, and per_effect_diff
#' @export
#'
#' @examples
check_pip_changes_per_effect <- function(PIP, PIP_prev, threshold = 1e-4, variant_names = NULL) {
  per_effect_diff <- PIP_prev - PIP
  idx <- which(rowSums(abs(per_effect_diff) > threshold) > 0)
  
  list(
    changed_idx      = idx,
    changed_variants = if (!is.null(variant_names)) variant_names[idx] else idx,
    per_effect_diff  = per_effect_diff
  )
}