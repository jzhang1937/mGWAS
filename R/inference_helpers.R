jaccard_sim <- function(mat) {
  mat <- mat > 0
  inter <- tcrossprod(mat)  
  sums <- rowSums(mat)
  union <- outer(sums, sums, "+") - inter
  inter / union  
}


