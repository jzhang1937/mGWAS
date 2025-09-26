#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
oracle_treeWAS_wrapper <- function(data, method) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  snps.rec <- data$X.rec
  phen.rec <- data$y.rec
  results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                               tree = tree,
                                               test = method,
                                               snps.reconstruction = snps.rec,
                                               phen.reconstruction = phen.rec,
                                               plot.tree = FALSE,
                                               plot.manhattan = FALSE,
                                               plot.null.dist = FALSE,
                                               plot.null.dist.pairs = FALSE))
  results$dat <- NULL
  results$nonnulls <- data$snps.assoc
  results
}

#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
treeWAS_wrapper <- function(data, method) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                               tree = tree,
                                               test = method,
                                               plot.tree = FALSE,
                                               plot.manhattan = FALSE,
                                               plot.null.dist = FALSE,
                                               plot.null.dist.pairs = FALSE))
  results$dat <- NULL
  results$nonnulls <- data$snps.assoc
  results
}

pyseer_wrapper <- function(data, method, tmpdir = "pyseer_temp",  
                           output = "pyseer_out.txt",
                           extra_args = NULL,
                           keep_files = FALSE) {
  X <- data$X
  y <- data$y
  n <- length(y)
  tree <- data$tree
  method <- match.arg(method)
  
  # make sure temp directory exists
  if (!dir.exists(tmpdir)) {
    dir.create(tmpdir, recursive = TRUE)
  }
  
  # sample names check 
  if (is.null(rownames(X))) rownames(X) <- paste0("sample_", seq(1:n))
  if (is.null(names(y))) rownames(X)
  if (!all(names(y) %in% rownames(X))) {
    stop("Sample names in y must match rownames of X")
  }
  
  samples <- names(y)
  
  # make phenotype file
  pheno_file <- file.path(tmpdir, "phenotypes.txt")
  pheno_df <- data.frame(sample = samples, phenotype = y[samples],
                         stringsAsFactors = FALSE)
  write.table(pheno_df, pheno_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  # make genotype file
  geno_file <- file.path(tmpdir, "variants.pres.tsv")
  mat <- t(X[samples, , drop = FALSE])   
  geno_out <- data.frame(variant = rownames(mat), mat, check.names = FALSE)
  write.table(geno_out, geno_file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  
  # distances / similarities from tree 
  dist_file <- sim_file <- NULL
  if (!is.null(tree)) {
    tr <- if (file.exists(tree)) ape::read.tree(tree) else ape::read.tree(text = tree)
    tr <- ape::keep.tip(tr, samples)
    
    if (method == "fixed") {
      D <- ape::cophenetic.phylo(tr)[samples, samples]
      dist_file <- file.path(tmpdir, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = TRUE)
    }
    
    if (method == "mixed") {
      D <- ape::cophenetic.phylo[samples, samples]
      sim <- exp(-D)   # similarity from distances
      sim_file <- file.path(tmpdir, "similarities.txt")
      write.table(sim, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = TRUE)
    }
  }
  
  # make pyseer command
  cmd <- c("pyseer",
           paste0("--phenotypes ", pheno_file),
           paste0("--pres ", geno_file),
           paste0("--output-pattern ", file.path(tmpdir, output)))
  
  if (!is.null(dist_file))  cmd <- c(cmd, paste0("--distances ", dist_file))
  if (!is.null(sim_file))   cmd <- c(cmd, paste0("--similarities ", sim_file), "--lmm")
  if (!is.null(extra_args)) cmd <- c(cmd, extra_args)
  
  cmd_str <- paste(cmd, collapse = " ")
  message("Running: ", cmd_str)
  system(cmd_str)
  
  # get results
  res <- read.table(file.path(tmpdir, output), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # remove files
  if (!keep_files) unlink(tmpdir, recursive = TRUE)
  
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
hogwash_wrapper <- function(data) {
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
  # run hogwash
  hogwash::hogwash(pheno = as.matrix(phen), 
          geno = as.matrix(snps),  
          tree = tree,      
          file_name = timestamp,
          dir = "hogwash",
          perm = 10000,
          fdr = 0.1,
          test = "both")
  if (binary) {
    synchronous <- read.rda(paste0("hogwash/hogwash_synchronous_", timestamp, ".rda"))
    convergence <- read.rda(paste0("hogwash/hogwash_convergence_", timestamp, ".rda"))
    unlink("hogwash/*", recursive = TRUE)
    return(list(synchronous = synchronous, convergence = convergence))
  } else {
    continuous <- read.rda(paste0("hogwash/hogwash_continuous_", timestamp, ".rda"))
    unlink("hogwash/*", recursive = TRUE)
    return(continuous)
  }
  
}

bugwas_wrapper <- function(data) {
  
}