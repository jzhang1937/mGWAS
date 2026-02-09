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