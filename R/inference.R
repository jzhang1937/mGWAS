#' Title
#'
#' @param data a list object that contains an X matrix, y vector, phylogenetic tree,
#' and a number of mutations vector n.mts.
#' @param method vector of scores to compute in treeWAS
#' @param n.snps.sim number of samples to use
#'
#' @returns output from treeWAS
#' @export
#'
#' @examples
oracle_treeWAS_wrapper <- function(data, method = c("terminal", "simultaneous", "subsequent"), n.snps.sim = 10000) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  n.mts <- data$n.mts
  n.subs <- as.vector(table(data$n.mts))
  results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                               tree = tree, n.subs = n.subs,
                                               test = method,
                                               n.snps.sim = n.snps.sim,
                                               plot.tree = FALSE,
                                               plot.manhattan = FALSE,
                                               plot.null.dist = FALSE,
                                               plot.null.dist.pairs = FALSE))
  results$dat <- NULL
  results
}

#' Title
#'
#' @param data a list object that contains an X matrix, y vector, phylogenetic tree.
#' The rownames of X and names of y should align, and be a subset of tree$tip.label
#' @param method vector of scores to compute in treeWAS
#' @param n.snps.sim number of samples to use
#'
#' @returns output from treeWAS
#' @export
#'
#' @examples
treeWAS_tree_wrapper <- function(data, method = c("terminal", "simultaneous", "subsequent"), n.snps.sim = 10000) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                               tree = tree,
                                               test = method,
                                               n.snps.sim = n.snps.sim,
                                               plot.tree = FALSE,
                                               plot.manhattan = FALSE,
                                               plot.null.dist = FALSE,
                                               plot.null.dist.pairs = FALSE))
  results$dat <- NULL
  results
}

#' Title
#'
#' @param data a list object that contains an X matrix, y vector, phylogenetic tree (optional). 
#' The rownames of X and names of y should align, and be a subset of tree$tip.label if tree is supplied.
#' @param method vector of scores to compute in treeWAS
#' @param n.snps.sim number of samples to use
#'
#' @returns list object with corrected treeWAS p-values
#' @export
#'
#' @examples
treeWAS_wrapper <- function(data, method = c("terminal", "simultaneous", "subsequent"), n.snps.sim = 10000) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  if (is.null(tree)) {
    results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                                 test = method,
                                                 n.snps.sim = n.snps.sim,
                                                 plot.tree = FALSE,
                                                 plot.manhattan = FALSE,
                                                 plot.null.dist = FALSE,
                                                 plot.null.dist.pairs = FALSE))
  } else {
    tree_pruned <- tree_prune(tree, names(phen))
    results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                                 tree = tree,
                                                 test = method,
                                                 n.snps.sim = n.snps.sim,
                                                 plot.tree = FALSE,
                                                 plot.manhattan = FALSE,
                                                 plot.null.dist = FALSE,
                                                 plot.null.dist.pairs = FALSE))
  }
  results$dat <- NULL
  process_treeWAS(results)
}

#' Title
#'
#' @param data a list object that contains an X matrix, y vector, phylogenetic tree (optional). 
#' The rownames of X and names of y should align, and be a subset of tree$tip.label if tree is supplied.
#' @param method "fixed" or "mixed", no default.
#' @param tmpdir a path to directory to put temporary files.
#' @param output name of output file.
#' @param conda_bin path to the conda executable
#' @param pyseer_env name of conda environment
#' @param extra_args 
#' @param keep_files whether to keep the temp files, default FALSE.
#' @param max_dimensions if using method = "fixed", how many MDS dimensions to use.
#' @param phylo whether to use tree similarity/distance (tree must be provided).
#' @param jaccard whether to use jaccard similarity/distance.
#' @param hamming whether to use hamming similarity/distance.
#' @param normalized whether to use a normalized XX^T similarity/distance.
#'
#' @returns pyseer results output
#' @export
#'
#' @examples
pyseer_wrapper <- function(data, method, tmpdir = "pyseer_temp",  
                           output = "pyseer_out.txt",
                           conda_bin,
                           pyseer_env,
                           max_dimensions = 10,
                           phylo = FALSE,
                           jaccard = FALSE,
                           hamming = FALSE,
                           normalized = FALSE,
                           extra_args = NULL,
                           keep_files = FALSE) {
  X <- data$X
  y <- data$y
  n <- length(y)
  tree <- data$tree
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  
  # make sure temp directory exists
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE)
  }
  tmpdir_time <- file.path(tmpdir, timestamp)
  dir.create(tmpdir_time, recursive = TRUE)
  
  # sample names standard
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  
  # make phenotype file
  pheno_file <- file.path(tmpdir_time, "phenotypes.txt")
  pheno_df <- data.frame(sample = samples, phenotype = y[samples],
                         stringsAsFactors = FALSE)
  write.table(pheno_df, pheno_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  # make genotype file
  geno_file <- file.path(tmpdir_time, "variants.pres.tsv")
  mat <- t(X[samples, , drop = FALSE])   
  geno_out <- data.frame(Gene = rownames(mat), mat, check.names = FALSE)
  write.table(geno_out, geno_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  # make sample list file
  sample_file <- file.path(tmpdir_time, "sample_list.txt")
  writeLines(samples, sample_file)
  
  # distances / similarities from tree 
  dist_file <- sim_file <- NULL
  if (!is.null(tree) & phylo) {
    tr <- ape::keep.tip(tree, samples)
    
    if (method == "fixed") {
      D <- ape::cophenetic.phylo(tr)[samples, samples]
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    } else if (method == "mixed") {
      sim <- ape::vcv.phylo(tr)[samples, samples]
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(sim, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    } 
  } else if (!phylo & jaccard) {
    jaccard_sim <- mGWAS:::jaccard_sim(X)
    if (method == "fixed") {
      # compute distance with jaccard
      D <- 1 - jaccard_sim
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
    else if (method == "mixed") {
      # compute similarity with jaccard
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(jaccard_sim, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    }
  } else if (!phylo & hamming) {
    G = 2 * X - 1
    K <- tcrossprod(G)
    K_norm <- K / ncol(G)
    if (method == "fixed") {
      # compute distance with hammming
      D <- 1 - K_norm
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
    else if (method == "mixed") {
      # compute similarity with hamming
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(K_norm, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    }
  } else if (!phylo & normalized) {
    # Allele frequencies
    freq <- colMeans(X, na.rm = TRUE) / 2
    # Center and scale
    Z <- sweep(X, 2, freq, "-")
    Z <- sweep(Z, 2, sqrt(freq * (1 - freq)), "/")
    
    # Drop degenerate columns
    keep <- is.finite(colSums(Z))
    Z <- Z[, keep, drop = FALSE]
    
    # GRM
    K <- tcrossprod(Z) / ncol(Z)
    
    if (method == "fixed") {
      # compute kinship with normalized 1 - GG^T
      D2 <- outer(diag(K), diag(K), "+") - 2 * K
      D  <- sqrt(pmax(D2, 0))   # numerical safety
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
    else if (method == "mixed") {
      # compute kinship
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(K, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
  } else {
    if (method == "fixed") {
      # compute kinship with normalized 1 - GG^T
      K <- tcrossprod(X)
      K_norm <- K / ncol(X)
      D <- 1 - K_norm
      # set diagonal to 0
      diag(D) <- 0
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
    else if (method == "mixed") {
      # compute kinship
      K <- X %*% t(X)
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(K, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
      
    }
  }
  
  # make pyseer command
  cmd <- c("pyseer",
           paste0("--phenotypes ", pheno_file),
           paste0("--pres ", geno_file))
  
  if (method == "fixed") {cmd <- c(cmd, paste0("--distances ", dist_file, " --max-dimensions ", max_dimensions)) }
  if (method == "mixed") {cmd <- c(cmd, paste0("--similarity ", sim_file), "--lmm")}
  if (method == "elasticnet") cmd <- c(cmd, "--wg enet --alpha 1")
  # if (!is.null(extra_args)) cmd <- c(cmd, extra_args)
  
  cmd_str <- paste(cmd, collapse = " ")
  # prepend conda run
  full_cmd <- paste(conda_bin, "run -n", pyseer_env, cmd_str,
                    paste0("> ", file.path(tmpdir_time, output)))
  
  message("Running: ", full_cmd)
  system(full_cmd)
  
  # get results
  res <- read.table(file.path(tmpdir_time, output), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # remove files
  if (!keep_files) unlink(file.path(tmpdir_time), recursive = TRUE)
  
  return(res)
}

#' Title
#'
#' @param data 
#' @param method 
#' @param tmpdir 
#' @param output 
#' @param conda_bin 
#' @param pyseer_env 
#' @param extra_args 
#' @param keep_files 
#'
#' @returns
#' @export
#'
#' @examples
scoary_wrapper <- function(data, tmpdir = "scoary_temp",  
                           N = 10000,
                           conda_bin,
                           scoary_env,
                           threads = 1,
                           low.var = NULL,
                           extra_args = NULL,
                           keep_files = FALSE) {
  X <- data$X
  y <- data$y
  
  n <- length(y)
  tree <- data$tree
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  
  # make sure temp directory exists
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE)
  }
  tmpdir_time <- file.path(tmpdir, timestamp)
  dir.create(tmpdir_time, recursive = TRUE)
  
  # sample names check 
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match colnames of X")
  }
  
  samples <- names(y)
  
  # make phenotype file
  pheno_file <- file.path(tmpdir_time, "pheno.csv")
  pheno_df <- data.frame(Name = samples, phenotype = y[samples],
                         stringsAsFactors = FALSE)
  write.csv(pheno_df, pheno_file, row.names = FALSE, quote = FALSE)
  
  # make genotype file
  geno_file <- file.path(tmpdir_time, "geno.csv")
  mat <- t(X[samples, , drop = FALSE])   
  geno_out <- data.frame(sample = rownames(mat), mat, check.names = FALSE)
  write.csv(geno_out, geno_file, row.names = FALSE, quote = FALSE)
  
  # tree if tree is nonnull
  if (!is.null(tree)) {
    # Write tree file
    tree_file <- file.path(tmpdir_time, "tree.nwk")
    ape::write.tree(tree, file = tree_file)
   
  }
  
  # make pyseer command
  cmd <- c("scoary",
           paste0("-t ", pheno_file),
           paste0("-g ", geno_file),
           paste0("-o ", tmpdir_time),
           paste0("-p ", 1.0),
           paste0("-c I B BH P"),
           paste0("-e ", N),
           paste0("-s ", 2),
           paste0("--threads ", threads))
  
  if (!is.null(tree_file))  cmd <- c(cmd, paste0("-n ", tree_file))
  
  cmd_str <- paste(cmd, collapse = " ")
  # prepend conda run
  full_cmd <- paste(conda_bin, "run -n", scoary_env, cmd_str)
  
  message("Running: ", full_cmd)
  system(full_cmd)
  
  # get results
  # list all results CSVs
  results_files <- list.files(tmpdir_time, pattern = "\\.results\\.csv$", full.names = TRUE)
  
  # trait name
  trait_name <- "phenotype"
  
  # match only the filename (without directory) with the trait
  trait_file <- results_files[grepl(paste0("^", trait_name), basename(results_files))]
  
  res <- read.csv(trait_file[1], stringsAsFactors = FALSE)
  
  # remove files
  if (!keep_files) unlink(file.path(tmpdir_time), recursive = TRUE)
  
  return(res)
}

#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
hogwash_wrapper <- function(data, method = "synchronous") {
  snps <- data$X
  phen <- data$y
  binary <- all(y %in% c(0, 1))
  tree <- data$tree
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # match names
  snps <- snps[tree$tip.label, , drop = FALSE]
  phen <- phen[tree$tip.label]
  # create directory if it doesn't exist
  if (!dir.exists("hogwash")) {
    dir.create("hogwash", recursive = TRUE)
  }
  # make labels
  tree$node.label <- rep(100, tree$Nnode)
  
  if (method == "synchronous") {
    # run hogwash
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "synchronous")
    load(paste0("hogwash/hogwash_synchronous_", timestamp, ".rda"))
    synchronous <- hogwash_synchronous
    return(synchronous)
  } else if (method == "phyc") {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "phyc")
    load(paste0("hogwash/hogwash_phyc_", timestamp, ".rda"))
    phyc <- hogwash_phyc
    unlink(paste0("hogwash/hogwash_phyc_", timestamp, "*"), recursive = TRUE)
    return(phyc)
  } else {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1)
    load(paste0("hogwash/hogwash_continuous_", timestamp, ".rda"))
    continuous <- hogwash_continuous
    unlink(paste0("hogwash/hogwash_continuous_", timestamp, "*"), recursive = TRUE)
    return(continuous)
  }
  
}

#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
one_match <- function(data, dist = "phylo", D = NA, cutoff = 1) {
  snps <- data$X
  phen <- data$y
  binary <- all(y %in% c(0, 1))
  tree <- data$tree
  if (is.na(D)) {
    if (dist == "phylo") {
      
    } else if (dist == "hamming") {
      
    } else if (dist == "jaccard") {
      
    }
  }
  
  
}

#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
minor_match <- function(data, method = "synchronous") {
  snps <- data$X
  phen <- data$y
  binary <- all(y %in% c(0, 1))
  tree <- data$tree
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # match names
  snps <- snps[tree$tip.label, , drop = FALSE]
  phen <- phen[tree$tip.label]
  # create directory if it doesn't exist
  if (!dir.exists("hogwash")) {
    dir.create("hogwash", recursive = TRUE)
  }
  # make labels
  tree$node.label <- rep(100, tree$Nnode)
  
  if (method == "synchronous") {
    # run hogwash
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "synchronous")
    load(paste0("hogwash/hogwash_synchronous_", timestamp, ".rda"))
    synchronous <- hogwash_synchronous
    return(synchronous)
  } else if (method == "phyc") {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "phyc")
    load(paste0("hogwash/hogwash_phyc_", timestamp, ".rda"))
    phyc <- hogwash_phyc
    unlink(paste0("hogwash/hogwash_phyc_", timestamp, "*"), recursive = TRUE)
    return(phyc)
  } else {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1)
    load(paste0("hogwash/hogwash_continuous_", timestamp, ".rda"))
    continuous <- hogwash_continuous
    unlink(paste0("hogwash/hogwash_continuous_", timestamp, "*"), recursive = TRUE)
    return(continuous)
  }
  
}

#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
full_match <- function(data, method = "synchronous") {
  snps <- data$X
  phen <- data$y
  binary <- all(y %in% c(0, 1))
  tree <- data$tree
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  # match names
  snps <- snps[tree$tip.label, , drop = FALSE]
  phen <- phen[tree$tip.label]
  # create directory if it doesn't exist
  if (!dir.exists("hogwash")) {
    dir.create("hogwash", recursive = TRUE)
  }
  # make labels
  tree$node.label <- rep(100, tree$Nnode)
  
  if (method == "synchronous") {
    # run hogwash
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "synchronous")
    load(paste0("hogwash/hogwash_synchronous_", timestamp, ".rda"))
    synchronous <- hogwash_synchronous
    return(synchronous)
  } else if (method == "phyc") {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1,
                     test = "phyc")
    load(paste0("hogwash/hogwash_phyc_", timestamp, ".rda"))
    phyc <- hogwash_phyc
    unlink(paste0("hogwash/hogwash_phyc_", timestamp, "*"), recursive = TRUE)
    return(phyc)
  } else {
    hogwash::hogwash(pheno = as.matrix(phen), 
                     geno = as.matrix(snps),  
                     tree = tree,      
                     file_name = timestamp,
                     dir = "hogwash",
                     perm = 10000,
                     bootstrap = 0.3,
                     fdr = 0.1)
    load(paste0("hogwash/hogwash_continuous_", timestamp, ".rda"))
    continuous <- hogwash_continuous
    unlink(paste0("hogwash/hogwash_continuous_", timestamp, "*"), recursive = TRUE)
    return(continuous)
  }
  
}

#' Title 
#'
#' @param data a list object that contains an X matrix, y vector, phylogenetic tree (optional). 
#' The rownames of X and names of y should align, and be a subset of tree$tip.label if tree is supplied.
#' @param phylo whether to use tree similarity/distance (tree must be provided).
#' @param jaccard whether to use jaccard similarity/distance.
#' @param hamming whether to use hamming similarity/distance.
#' @param normalized whether to use a normalized XX^T similarity/distance.
#' @param tmpdir path to temp directory to store output
#' @param output name of output file.
#' @param covar covariate matrix.
#'
#' @returns GMMAT results output
#' @export
#'
#' @examples
gmmat_binary <- function(data, covar = NULL, 
                             phylo = FALSE,
                             jaccard = FALSE,
                             hamming = FALSE,
                             normalized = FALSE,
                             tmpdir = "gmmat_temp",  
                             output = "gmmat_out.txt") {
  X <- data$X
  y <- data$y
  tree <- data$tree
  
  stopifnot(length(y) == nrow(X))
  
  n <- length(y)
  p <- ncol(X)
  
  # sample names standard
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
  
  # make sure temp directory exists
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE)
  }
  tmpdir_time <- file.path(tmpdir, timestamp)
  dir.create(tmpdir_time, recursive = TRUE)
  
  # Make data frame
  dat <- data.frame(
    id = seq_len(n),
    y = as.numeric(y)
  )
  
  if (!is.null(covar)) {
    covar <- as.data.frame(covar)
    dat <- cbind(dat, covar)
  }
  
  # Compute similarity K
  if (!is.null(tree) & phylo) {
    tr <- ape::keep.tip(tree, samples)
    K <- ape::vcv.phylo(tr)[samples, samples]
  } else if (!phylo & jaccard) {
    K <- mGWAS:::jaccard_sim(X)
  } else if (!phylo & hamming) {
    G = 2 * X - 1
    K_raw <- tcrossprod(G)
    K <- K_raw / ncol(G)
  } else if (!phylo & normalized) {
    # Allele frequencies
    freq <- colMeans(X, na.rm = TRUE) / 2
    # Center and scale
    Z <- sweep(X, 2, freq, "-")
    Z <- sweep(Z, 2, sqrt(freq * (1 - freq)), "/")
    
    # Drop degenerate columns
    keep <- is.finite(colSums(Z))
    Z <- Z[, keep, drop = FALSE]
    
    # GRM
    K <- tcrossprod(Z) / ncol(Z)
    
  } else {
    K <- tcrossprod(X)
  }
  
  rownames(K) <- colnames(K) <- dat$id
  
  # make genotype file
  geno_file <- file.path(tmpdir_time, "genotypes.txt")
  
  # Transpose so rows = variants, columns = samples
  mat <- t(X[samples, , drop = FALSE])  
  
  # Make a data.frame, optional to add a variant identifier column
  geno_out <- data.frame(
    variant = rownames(mat),  # optional, can be first column to label variants
    mat,
    check.names = FALSE
  )
  
  write.table(
    geno_out,
    file = geno_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  # Fit null logistic mixed model
  fixed <- if (is.null(covar)) {
    y ~ 1
  } else {
    reformulate(colnames(covar), response = "y")
  }
  
  message("Fitting null model...")
  nullmod <- GMMAT::glmmkin(
    fixed = fixed,
    data = dat,
    kins = list(K),
    id = "id",
    family = binomial(link = "logit")
  )
  
  # Score test for each variant
  message("Running score tests...")
  res <- GMMAT::glmm.score(
    nullmod,
    infile = geno_file,
    outfile = file.path(tmpdir_time, output)
  )
  
  results <- read.table(file.path(tmpdir_time, output), header = TRUE)
  results <- results[-1, ]
  
  # remove files
  unlink(file.path(tmpdir_time), recursive = TRUE)
  return(results)
}

#' GCM doubly robust variable selection function
#'
#' @param data A list with an n x p data matrix X and n x 1 vector y
#' @param fit_X A function that inputs X and outputs a function
#' that inputs X, j, B and outputs B conditional resamples of
#' variable j.
#' @param fit_y_given_X A function that inputs X and y and family and
#' outputs a function that inputs X and outputs the fitted
#' E(y|X)
#' @param family A string specifying the exponential family, e.g. "gaussian", "bernoulli"
#' @param K The number of folds for crossfitting
#' @param B The number of resamples
#' @param linear whether the fit_y_given_X is linear in X
#' @param multiple_correction How to adjust the p values for multiplicity.
#' @param seed Seed to set for data split and fits of Y|X and X
#' @param alpha Desired FDP level.
#'
#' @return A list of p p-values, one for each variable
#' @export
tower_gcm <- function(data, fit_X, fit_y_given_X, 
                     K = 5, alpha = 0.1,
                     multiple_correction = "BH", seed = 1234){
  # Extract the data.
  X <- data$X
  y <- data$y
  
  # number of samples and predictors
  n <- nrow(X)
  p <- ncol(X)
  
  # Create K roughly equally size folds
  folds <- sample(x = 1:K, size = n, replace = TRUE) |>
    R.utils::withSeed(seed = seed)
  
  # Matrix of product of residuals across folds
  product_residuals <- matrix(0,nrow = n, ncol = p)
  
  
  #Perform K-fold cross-fitting
  for(k in 1:K){
    print(paste0("Fitting Nuisance Fold ", k))
    # Segment data by fold using the which() function
    evalIndexes <- which(folds==k,arr.ind=TRUE)
    # evaluation data
    evalDataX <- X[evalIndexes, ]
    evalDataY <- y[evalIndexes]
    # nuisance training data
    if (K > 1) {
      nuisDataX <- X[-evalIndexes, ]
      nuisDataY <- y[-evalIndexes]
    } else if (K == 1) {
      nuisDataX <- evalDataX
      nuisDataY <- evalDataY
    }
    
    # Concatenate all of the training data
    data_nuis <- list(X = nuisDataX, y = nuisDataY)
    
    # number of evaluation samples
    nk <- length(evalDataY)
    
    # y given X fit
    y_given_X_hat <- fit_y_given_X(data_nuis) |>
      R.utils::withSeed(seed = seed)
    
    # X fit
    X_hat <- fit_X(data_nuis) |>
      R.utils::withSeed(seed = seed)
    
    # print update
    print(paste0("Fitting Nuisance Complete, Fold: ", k))
    
    # Iterate over each predictor variable
    cat(sprintf("Running association test for each variable...\n"))
    # Iterate over each predictor variable
    for (j in 1:p) {

      # compute E[eta^hat_j(evalDataX)|X_-j]
      cond_mean_xj <- X_hat$conditional_mean(j, evalDataX)
      # round to be between 0 and 1
      cond_probs_xj <- pmin(pmax(cond_mean_xj, 0), 1)
      condProbs <- cbind(as.vector(1 - cond_probs_xj), as.vector(cond_probs_xj))
      
      evalDataX0j <- evalDataX
      evalDataX0j[,j] <- 0
      evalDataX1j <- evalDataX
      evalDataX1j[,j] <- 1
      
      # 3) Compute f(X) for each slice (dimension 3)
      possible_y_X <- cbind(y_given_X_hat(evalDataX0j), y_given_X_hat(evalDataX1j))
      
      # Compute the conditional mean (length n vector)
      cond_mean_yj <- rowSums(possible_y_X * condProbs)  # (n x num_obs) * (n x num_obs) → sum along num_obs → (n)
      
      # compute the contribution to the test statistic for all subjects from
      # the current evaluation fold. (X_j - E(X_j|X_-j))(Y - E(Y - E(Y|X_-j)))
      product_residuals[evalIndexes,j] <- (evalDataX[,j] - cond_mean_xj) *
        (evalDataY - cond_mean_yj)
      
    }
    
  }
  # standard deviation estimator
  sd_hat = apply(X = product_residuals, MARGIN = 2, FUN = stats::sd)
  
  # calculate test statistics and p-values.
  test_stats <- colSums(product_residuals)/(sqrt(n)*sd_hat)
  test_stats[is.na(test_stats)] <- 0
  test_stats[abs(colMeans(product_residuals)) < .Machine$double.eps] <- 0
  p_values <- 2*(1-stats::pnorm(abs(test_stats))) # two-sided p-value
  #stats::pnorm(test_stats, lower.tail = FALSE)  one-sided p-value
  
  cat(sprintf("Done.\n"))
  
  # multiple testing correction
  p_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (multiple_correction %in% p_methods) {
    p_values_correction <- stats::p.adjust(p_values, multiple_correction)
    selected <- unname(which(p_values_correction < alpha))
  } else {
    selected = which(p_values < alpha)
  } # will add more options later.
  
  # return
  return(list(p_values = p_values, nonnulls = selected))
  
}

#' Title
#'
#' @param data 
#' @param fit_X 
#' @param fit_y_given_X 
#' @param family 
#'
#' @returns
#' @export
#'
#' @examples
spaVS <- function(data, fit_X, fit_y_given_X, family = "binomial") {
  
}



