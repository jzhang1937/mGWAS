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
