generate_data_treeWAS <- function(n.ind, n.snps, n.snps.assoc, assoc.prob) {
  data <- treeWAS::coalescent.sim(n.ind = n.ind, n.snps = n.snps,
                                  n.snps.assoc = n.snps.assoc,
                                  assoc.prob = assoc.prob)
  return(data)
}

generate_data_treeWAS_mod <- function(n.ind, n.snps, n.snps.assoc, assoc.prob, ground_truth) {
  data <- coalescent.sim.mod(n.ind = n.ind, n.snps = n.snps,
                                  n.snps.assoc = n.snps.assoc,
                                  assoc.prob = assoc.prob,
                             ground.truth = ground_truth)
  return(data)
}

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
