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




#' Title
#'
#' @param data a list containing a vector of outcomes y and data matrix X.
#' @param L number of modeled causal effects
#' @param V precomputed p x p matrix of eigenvectors of X'X
#' @param Dsq precomputed length-p vector of eigenvalues of X'X
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param method one of {'moments','MLE'} (sigma^2,tau^2) are estimated using method-of-moments or MLE
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns List containing PIP -- p x L matrix of PIPs, individually for each effect mu -- p x L matrix of posterior means conditional on causal omega -- p x L matrix of posterior precisions conditional on causal lbf -- length-L array of log-Bayes-factor for each effect lbf_variable -- p x L matrix of per-variable log-Bayes-factors ssq -- length-L array of final effect size variances s^2 sigmasq -- final value of sigma^2 tausq -- final value of tau^2 alpha -- length-p array of posterior means of infinitesimal effects
#' @export
#'
#' @examples
susieinf <- function(data, L, V = NULL, Dsq = NULL,
                    est_ssq = TRUE, ssq = NULL, ssq_range = c(0,1), pi0 = NULL,
                    est_sigmasq = TRUE, est_tausq = TRUE,
                    sigmasq = 1, tausq = 0,
                    method = "moments",
                    sigmasq_range = NULL, tausq_range = NULL,
                    PIP = NULL, mu = NULL,
                    maxiter = 100, PIP_tol = 1e-3, verbose = TRUE) {
  y <- data$y
  y <- y - mean(y)
  X <- data$X
  # center columns
  X <- sweep(X, 2, colMeans(X), "-")
  col_sd <- sqrt(colMeans(X^2))   
  X <- sweep(X, 2, col_sd, "/")
  n <- nrow(X); p <- ncol(X)
  
  # sample names standard
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  variants <- colnames(X)
  
  # LD matrix
  LD <- crossprod(X) / n
  # z <- t(X) %*% y / sqrt(n)
  z <- as.vector(crossprod(X, y) / sqrt(n))
  meansq <- sum(y^2) / n
  
  if ((is.null(V) || is.null(Dsq)) && is.null(LD))
    stop("Missing LD")
  
  if (is.null(V) || is.null(Dsq)) {
    eig <- eigen(LD, symmetric = TRUE)
    V <- eig$vectors
    Dsq <- pmax(n * eig$values, 0)
  } else {
    Dsq <- pmax(Dsq, 0)
  }
  
  Xty <- sqrt(n) * z
  VtXty <- crossprod(V, Xty)
  yty <- n * meansq
  
  var <- tausq * Dsq + sigmasq
  # diagXtOmegaX <- colSums(V^2 * (Dsq / var))
  diagXtOmegaX <- colSums(t(V^2) * (Dsq / var))
  XtOmegay <- as.vector(V %*% (VtXty / var))
  
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, p, L)
  if (is.null(mu)) mu <- matrix(0, p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf <- numeric(L)
  omega <- matrix(0, nrow = p, ncol = L)
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  # omega <- diagXtOmegaX + 1 / ssq
  
  logpi0 <- if (is.null(pi0)) rep(log(1/p), p) else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  for (it in seq_len(maxiter)) {
    
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    for (l in seq_len(L)) {
      
      b <- rowSums(mu * PIP) - mu[,l] * PIP[,l]
      
      XtOmegaXb <- V %*% ((crossprod(V, b) * Dsq) / var)
      XtOmegar <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) +
            logpi0
          -logsumexp(val)
        }
        
        opt <- optimize(f, ssq_range, tol = 1e-5)
        
        if (is.finite(opt$objective) &&
            opt$minimum > ssq_range[1] &&
            opt$minimum < ssq_range[2]) {
          ssq[l] <- opt$minimum
        } else if (verbose) {
          cat(sprintf("WARNING: s^2 update failed at effect %d\n", l))
        }
      }
      
      omega[,l] <- diagXtOmegaX + 1 / ssq[l]
      mu[,l] <- XtOmegar / omega[,l]
      
      lbf_variable[,l] <- XtOmegar^2 / (2 * omega[,l]) -
        0.5 * log(omega[,l] * ssq[l])
      
      logPIP <- lbf_variable[,l] + logpi0
      lbf[l] <- logsumexp(logPIP)
      PIP[,l] <- exp(logPIP - lbf[l])
    }
    
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        res <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq,
                   VtXty, Xty, yty, est_sigmasq, est_tausq, verbose)
      } else {
        res <- MLE(PIP, mu, omega, sigmasq, tausq, n, V, Dsq,
                   VtXty, yty, est_sigmasq, est_tausq,
                   sigmasq_range, tausq_range, it, verbose)
      }
      sigmasq <- res$sigmasq; tausq <- res$tausq
      
      var <- tausq * Dsq + sigmasq
      # diagXtOmegaX <- colSums(V^2 * (Dsq / var))
      diagXtOmegaX <- colSums(t(V^2) * (Dsq / var))
      XtOmegay <- V %*% (VtXty / var)
    }
    
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("Max PIP diff:", PIP_diff, "\n")
    
    if (PIP_diff < PIP_tol) {
      converged <- TRUE
      break
    }
  }
  
  if (!exists("converged")) converged <- FALSE
  
  b <- rowSums(mu * PIP)
  XtOmegaXb <- V %*% ((crossprod(V, b) * Dsq) / var)
  XtOmegar <- XtOmegay - XtOmegaXb
  alpha <- tausq * XtOmegar
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(PIP = PIP, marginalPIP = marginalPIP, mu = mu, omega = omega,
       lbf = lbf, lbf_variable = lbf_variable,
       ssq = ssq, sigmasq = sigmasq,
       tausq = tausq, alpha = alpha,
       converged = converged, variants = variants)
}

#' Title
#'
#' @param data a list containing a vector of outcomes y, data matrix X, and optionally a tree
#' @param L number of modeled causal effects
#' @param K a similarity matrix
#' @param phylo whether to use the given tree to construct the similarity matrix K
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param method one of {'moments','MLE'} (sigma^2,tau^2) are estimated using method-of-moments or MLE
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns List containing PIP -- p x L matrix of PIPs, individually for each effect mu -- p x L matrix of posterior means conditional on causal omega -- p x L matrix of posterior precisions conditional on causal lbf -- length-L array of log-Bayes-factor for each effect lbf_variable -- p x L matrix of per-variable log-Bayes-factors ssq -- length-L array of final effect size variances s^2 sigmasq -- final value of sigma^2 tausq -- final value of tau^2 alpha -- length-p array of posterior means of infinitesimal effects
#' @export
#'
#' @examples
susieK <- function(data, L, K = NULL, phylo = TRUE,
                     est_ssq = TRUE, ssq = NULL, ssq_range = c(0,1), pi0 = NULL,
                     est_sigmasq = TRUE, est_tausq = TRUE,
                     sigmasq = 1, tausq = 0,
                     method = "moments",
                     sigmasq_range = NULL, tausq_range = NULL,
                     PIP = NULL, mu = NULL,
                     maxiter = 100, PIP_tol = 1e-3, verbose = TRUE) {
  y <- data$y
  y <- y - mean(y)
  X <- data$X
  # center columns
  X <- sweep(X, 2, colMeans(X), FUN = "-")
  # scale columns 
  X <- sweep(X, 2, apply(X, 2, sd), FUN = "/")
  tree <- data$tree
  n <- nrow(X); p <- ncol(X)
  
  # sample names standard
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  variants <- colnames(X)
  
  if (is.null(K) & !phylo) {
    K = X %*% t(X)
    message("K required if not using tree; using XX^T")
  }
  if (phylo) K <- ape::vcv.phylo(tree)[samples, samples]
  
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, p, L)
  if (is.null(mu)) mu <- matrix(0, p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf <- numeric(L)
  
  logpi0 <- if (is.null(pi0)) rep(log(1/p), p) else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  eig <- eigen(K, symmetric = TRUE)
  eigvals <- eig$values
  Q <- eig$vectors
  
  Z <- crossprod(Q, X)          # Q^T X
  y_tilde <- crossprod(Q, y)    # Q^T y
  
  w <- 1 / (tausq * eigvals + sigmasq)
  
  diagXtOmegaX <- colSums(Z^2 * w)
  XtOmegay <- crossprod(Z, w * y_tilde)
  
  # precompute invariants
  yty <- crossprod(y,y)
  Xty <- crossprod(X, y)
  
  trK <- sum(eigvals)
  trK2 <- sum(eigvals^2)
  
  yKy <- sum(eigvals * y_tilde^2)
  XTKy <- crossprod(Z, eigvals * y_tilde)
  
  XtX_diag <- colSums(Z^2)
  XTKX_diag <- colSums(Z^2 * eigvals)
  
  omega <- matrix(0, nrow = p, ncol = L)
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  # omega <- diagXtOmegaX + 1 / ssq
  
  for (it in seq_len(maxiter)) {
    
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    for (l in seq_len(L)) {
      
      b <- rowSums(mu * PIP) - mu[,l] * PIP[,l]
      
      XtOmegaXb <- crossprod(Z, w * (Z %*% b))
      XtOmegar <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) + logpi0
          -logsumexp(val)
        }
        opt <-  optimize(f, ssq_range, tol = 1e-5)
        ssq[l] <- opt$minimum
      }
      
      omega[,l] <- diagXtOmegaX + 1 / ssq[l]
      mu[,l] <- XtOmegar / omega[,l]
      
      lbf_variable[,l] <- XtOmegar^2 / (2 * omega[,l]) -
        0.5 * log(omega[,l] * ssq[l])
      
      logPIP <- lbf_variable[,l] + logpi0
      lbf[l] <- logsumexp(logPIP)
      PIP[,l] <- exp(logPIP - lbf[l])
    }
    
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        res <- MoMK(PIP, mu, omega, sigmasq, tausq, n, eigvals,
                    XtX_diag, XTKX_diag, Xty, XTKy,
                    yty, yKy, trK, trK2,
                    est_sigmasq, est_tausq, verbose)
        sigmasq <- res$sigmasq
        tausq <- res$tausq
      } else {
        stop("MLE not implemented for K")
      }
      
      w <- 1 / (tausq * eigvals + sigmasq)
      diagXtOmegaX <- colSums(Z^2 * w)
      XtOmegay <- crossprod(Z, w * y_tilde)
    }
    
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("Max PIP diff:", PIP_diff, "\n")
    
    if (PIP_diff < PIP_tol) {
      converged <- TRUE
      break
    }
  }
  
  if (!exists("converged")) converged <- FALSE
  
  b <- rowSums(mu * PIP)
  r <- y - X %*% b
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(PIP = PIP, marginalPIP = marginalPIP, mu = mu, omega = omega,
       lbf = lbf, lbf_variable = lbf_variable,
       ssq = ssq, sigmasq = sigmasq,
       tausq = tausq, converged = converged, variants = variants)
}

#' Title
#'
#' @param data a list containing a vector of outcomes y and data matrix X.
#' @param L number of modeled causal effects
#' @param V precomputed p x p matrix of eigenvectors of X'X
#' @param Dsq precomputed length-p vector of eigenvalues of X'X
#' @param est_ssq estimate prior effect size variances s^2 using MLE
#' @param ssq length-L initialization s^2 for each effect Default: 0.2 for every effect
#' @param ssq_range lower and upper bounds for each s^2, if estimated
#' @param pi0 length-p vector of prior causal probability for each SNP; must sum to 1 Default: 1/p for every SNP
#' @param est_sigmasq estimate variance sigma^2
#' @param est_tausq estimate both variances sigma^2 and tau^2
#' @param sigmasq initial value for sigma^2
#' @param tausq initial value for tau^2
#' @param method one of {'moments','MLE'} (sigma^2,tau^2) are estimated using method-of-moments or MLE
#' @param sigmasq_range lower and upper bounds for sigma^2, if estimated using MLE 
#' @param tausq_range lower and upper bounds for tau^2
#' @param PIP p x L initializations of PIPs Default: 1/#SNPs for each SNP and effect
#' @param mu p x L initializations of mu; Default: 0 for each SNP and effect
#' @param maxiter maximum number of SuSiE iterations
#' @param PIP_tol convergence threshold for PIP difference between iterations
#' @param verbose 
#'
#' @returns List containing PIP -- p x L matrix of PIPs, individually for each effect mu -- p x L matrix of posterior means conditional on causal omega -- p x L matrix of posterior precisions conditional on causal lbf -- length-L array of log-Bayes-factor for each effect lbf_variable -- p x L matrix of per-variable log-Bayes-factors ssq -- length-L array of final effect size variances s^2 sigmasq -- final value of sigma^2 tausq -- final value of tau^2 alpha -- length-p array of posterior means of infinitesimal effects
#' @export
#'
#' @examples
susieinf_wrapper <- function(data, L, V = NULL, Dsq = NULL,
                     est_ssq = TRUE, ssq = NULL, ssq_range = c(0,1), pi0 = NULL,
                     est_sigmasq = TRUE, est_tausq = TRUE,
                     sigmasq = 1, tausq = 0,
                     method = "moments",
                     sigmasq_range = NULL, tausq_range = NULL,
                     PIP = NULL, mu = NULL,
                     maxiter = 100, PIP_tol = 1e-3, verbose = TRUE) {
  y <- data$y
  y <- y - mean(y)
  X <- data$X
  # center columns
  X <- sweep(X, 2, colMeans(X), "-")
  col_sd <- sqrt(colMeans(X^2))   
  X <- sweep(X, 2, col_sd, "/")
  n <- nrow(X); p <- ncol(X)
  
  # sample names standard
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) names(y) <- rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  variants <- colnames(X)
  
  # LD matrix
  LD <- crossprod(X) / n
  # z <- t(X) %*% y / sqrt(n)
  z <- as.vector(crossprod(X, y) / sqrt(n))
  meansq <- sum(y^2) / n
  
  if ((is.null(V) || is.null(Dsq)) && is.null(LD))
    stop("Missing LD")
  
  if (is.null(V) || is.null(Dsq)) {
    eig <- eigen(LD, symmetric = TRUE)
    V <- eig$vectors
    Dsq <- pmax(n * eig$values, 0)
  } else {
    Dsq <- pmax(Dsq, 0)
  }
  
  Xty <- sqrt(n) * z
  VtXty <- crossprod(V, Xty)
  yty <- n * meansq
  
  var <- tausq * Dsq + sigmasq
  # diagXtOmegaX <- colSums(V^2 * (Dsq / var))
  diagXtOmegaX <- colSums(t(V^2) * (Dsq / var))
  XtOmegay <- as.vector(V %*% (VtXty / var))
  
  if (is.null(ssq)) ssq <- rep(0.2, L)
  if (is.null(PIP)) PIP <- matrix(1/p, p, L)
  if (is.null(mu)) mu <- matrix(0, p, L)
  
  lbf_variable <- matrix(0, p, L)
  lbf <- numeric(L)
  omega <- matrix(0, nrow = p, ncol = L)
  omega <- diagXtOmegaX + matrix(1 / ssq, nrow = p, ncol = L, byrow = TRUE)
  # omega <- diagXtOmegaX + 1 / ssq
  
  logpi0 <- if (is.null(pi0)) rep(log(1/p), p) else {
    out <- rep(-Inf, p)
    out[pi0 > 0] <- log(pi0[pi0 > 0])
    out
  }
  
  for (it in seq_len(maxiter)) {
    
    if (verbose) cat("Iteration", it, "\n")
    PIP_prev <- PIP
    
    for (l in seq_len(L)) {
      
      b <- rowSums(mu * PIP) - mu[,l] * PIP[,l]
      
      XtOmegaXb <- V %*% ((crossprod(V, b) * Dsq) / var)
      XtOmegar <- XtOmegay - XtOmegaXb
      
      if (est_ssq) {
        f <- function(x) {
          val <- -0.5 * log(1 + x * diagXtOmegaX) +
            x * XtOmegar^2 / (2 * (1 + x * diagXtOmegaX)) +
            logpi0
          -logsumexp(val)
        }
        
        opt <- optimize(f, ssq_range, tol = 1e-5)
        
        if (is.finite(opt$objective) &&
            opt$minimum > ssq_range[1] &&
            opt$minimum < ssq_range[2]) {
          ssq[l] <- opt$minimum
        } else if (verbose) {
          cat(sprintf("WARNING: s^2 update failed at effect %d\n", l))
        }
      }
      
      omega[,l] <- diagXtOmegaX + 1 / ssq[l]
      mu[,l] <- XtOmegar / omega[,l]
      
      lbf_variable[,l] <- XtOmegar^2 / (2 * omega[,l]) -
        0.5 * log(omega[,l] * ssq[l])
      
      logPIP <- lbf_variable[,l] + logpi0
      lbf[l] <- logsumexp(logPIP)
      PIP[,l] <- exp(logPIP - lbf[l])
    }
    
    if (est_sigmasq || est_tausq) {
      if (method == "moments") {
        res <- MoM(PIP, mu, omega, sigmasq, tausq, n, V, Dsq,
                   VtXty, Xty, yty, est_sigmasq, est_tausq, verbose)
      } else {
        res <- MLE(PIP, mu, omega, sigmasq, tausq, n, V, Dsq,
                   VtXty, yty, est_sigmasq, est_tausq,
                   sigmasq_range, tausq_range, it, verbose)
      }
      sigmasq <- res$sigmasq; tausq <- res$tausq
      
      var <- tausq * Dsq + sigmasq
      # diagXtOmegaX <- colSums(V^2 * (Dsq / var))
      diagXtOmegaX <- colSums(t(V^2) * (Dsq / var))
      XtOmegay <- V %*% (VtXty / var)
    }
    
    PIP_diff <- max(abs(PIP_prev - PIP))
    if (verbose) cat("Max PIP diff:", PIP_diff, "\n")
    
    if (PIP_diff < PIP_tol) {
      converged <- TRUE
      break
    }
  }
  
  if (!exists("converged")) converged <- FALSE
  
  b <- rowSums(mu * PIP)
  XtOmegaXb <- V %*% ((crossprod(V, b) * Dsq) / var)
  XtOmegar <- XtOmegay - XtOmegaXb
  alpha <- tausq * XtOmegar
  marginalPIP <- 1 - apply(1 - PIP, 1, prod)
  
  list(PIP = PIP, marginalPIP = marginalPIP, mu = mu, omega = omega,
       lbf = lbf, lbf_variable = lbf_variable,
       ssq = ssq, sigmasq = sigmasq,
       tausq = tausq, alpha = alpha,
       converged = converged, variants = variants)
}