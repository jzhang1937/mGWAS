#' Title
#'
#' @param n.ind 
#' @param n.snps 
#' @param n.snps.assoc 
#' @param assoc.prob 
#'
#' @returns
#' @export
#'
#' @examples
generate_data_treeWAS <- function(n.ind, n.snps, n.snps.assoc, assoc.prob) {
  data <- treeWAS::coalescent.sim(n.ind = n.ind, n.snps = n.snps,
                                  n.snps.assoc = n.snps.assoc,
                                  assoc.prob = assoc.prob)
  return(data)
}

#' Title
#'
#' @param n.ind 
#' @param n.snps 
#' @param n.snps.assoc 
#' @param assoc.prob 
#' @param set 
#' @param ground_truth 
#'
#' @returns
#' @export
#'
#' @examples
generate_data_treeWAS_mod <- function(n.ind, n.snps, n.snps.assoc, assoc.prob = 95, set = 1, ground_truth) {
  data <- coalescent.sim.mod(n.ind = n.ind, n.snps = n.snps,
                                  n.snps.assoc = n.snps.assoc,
                                  assoc.prob = assoc.prob, set = set,
                             ground.truth = ground_truth)
  return(data)
}
#' Title
#'
#' @param n 
#' @param p 
#' @param s 
#' @param joint_X 
#' @param y_given_X 
#' @param X_hyperparams 
#' @param y_given_X_hyperparams 
#' @param amplitude 
#' @param ground_truth 
#'
#' @returns
#' @export
#'
#' @examples
generate_data <- function(n, p, s, joint_X, y_given_X, X_hyperparams,
                          y_given_X_hyperparams, amplitude, ground_truth) {
  # get the list of arguments
  data_gen_args <- as.list(environment())
  
  # generate genome
  if (joint_X == "treeWAS") {
    assoc.prob <- as.integer(X_hyperparams$assoc.prob)
    set <- as.integer(X_hyperparams$set)
    sim.data <- generate_data_treeWAS_mod(n.ind = n, n.snps = p, n.snps.assoc = s,
                                          assoc.prob = assoc.prob, set = set,
                                          ground_truth = ground_truth)
    X <- sim.data$snps
    X.rec <- sim.data$snps.rec
    y <- sim.data$phen - 1
    y.rec <- sim.data$phen.rec
    tree <- sim.data$tree
    nonnulls <- ground_truth$nonnulls
    
  } else if (joint_X == "simurg") {
    ref <- 'ref_tutorial.fasta'
    prob.gene.gain <- as.numeric(X_hyperparams$prob.gene.gain)
    prob.gene.loss <- as.numeric(X_hyperparams$prob.gene.loss)
    sub.rate <- as.numeric(X_hyperparams$sub.rate)
    sim.data <- simurg::simpg(ref = ref, norg = n, ne = 1e10, C = p, u = prob.gene.gain, 
               v = prob.gene.loss, mu = sub.rate, write_by = 'genome', 
               dir_out = 'simurg_temp', force = TRUE, replace = TRUE)
    X <- sim.data$panmatrix
    X <- X[, apply(X, 2, var) > 0.02]
    p <- ncol(X)
    tree <- sim.data$coalescent
    tree$node.label <- rep(1, tree$Nnode)
    unlink("simurg_temp/*", recursive = TRUE)
    X.rec <- NULL
    nonnulls <- ground_truth$nonnulls
  } 
  # generate phenotype
  # set up the exponential family object
  family <- y_given_X_hyperparams$family
  family_object <- eval(parse(text = sprintf("stats::%s()", family)))
  
  # generate y's.
  if (y_given_X == "linear") {
    beta <- amplitude * (1:p %in% ground_truth$nonnulls)
    if (is.null(X.rec)) {
      X.rec <- X
    }
    y.rec <- katlabutils::generate_glm_response_data(
      X.rec,
      beta,
      family
    )
    y <- y.rec[1:n]
    names(y) <- rownames(X)
  } else if (y_given_X == "gam") {
    beta <- amplitude * (1:p %in% ground_truth$nonnulls)
    transforms <- y_given_X_hyperparams$transforms
    # Make transformed X
    if (is.null(X.rec)) {
      X.rec <- X
    }
    transformed_X <- t(apply(X.rec, MARGIN = 1, function(row) {
      sapply(1:ncol(X.rec), function(col) {
        transforms[[col]](row[col])
      })
    }))
    y <- katlabutils::generate_glm_response_data(
      transformed_X,
      beta,
      family
    )
    y <- y.rec[1:n]
    names(y) <- rownames(X)
  } else if (y_given_X == "interacted") {
    # nonnulls
    if (is.null(X.rec)) {
      X.rec <- X
    }
    X_subset <- X.rec[,ground_truth$nonnulls]
    X_interact <- X_subset
    # coefficients
    beta <- rep(1, ncol(X_subset))
    # highest level of interactions
    order <- y_given_X_hyperparams$order
    # proportion by which to decay signal for each level of interaction
    decay <- y_given_X_hyperparams$decay
    # transform to apply
    transform <- y_given_X_hyperparams$transform
    for (k in 2:order) {
      # Get all combinations of k columns
      combos <- combn(ncol(X_subset), k, simplify = FALSE)
      
      # Compute interaction terms for each combination
      interactions <- sapply(combos, function(combo) {
        apply(X_subset[, combo, drop = FALSE], 1, prod)
      })
      # update transformed X
      X_interact <- cbind(X_interact, interactions)
      # update coefficient vector
      beta <- c(beta, rep(1 * decay^(k-1), ncol(interactions)))
      
    }
    # apply transform
    transformed_X <- transform(X_interact %*% beta)
    # draw y
    y.rec <- katlabutils::generate_glm_response_data(
      transformed_X,
      amplitude,
      family
    )
    y <- y.rec[1:n]
    names(y) <- rownames(X)
  }
  
  data <- list(X = X, y = y, X.rec = X.rec, y.rec = y.rec, tree = tree, 
               nonnulls = nonnulls, data_gen_args = data_gen_args)
}

