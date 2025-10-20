#' Title
#'
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
oracle_treeWAS_wrapper <- function(data, method, n.snps.sim = 10000) {
  snps <- data$X
  phen <- data$y
  tree <- data$tree
  snps.rec <- data$X.rec
  phen.rec <- data$y.rec
  n.mts <- data$n.mts
  n.subs <- as.vector(table(sim.data$n.mts))
  results <- suppressWarnings(treeWAS::treeWAS(snps = snps, phen = phen,
                                               tree = tree, n.subs = n.subs,
                                               test = method,
                                               snps.reconstruction = snps.rec,
                                               phen.reconstruction = phen.rec,
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
#' @param data 
#'
#' @returns
#' @export
#'
#' @examples
treeWAS_wrapper <- function(data, method, n.snps.sim = 10000) {
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
pyseer_wrapper <- function(data, method, tmpdir = "pyseer_temp",  
                           output = "pyseer_out.txt",
                           conda_bin,
                           pyseer_env,
                           phylo = FALSE,
                           pyseer_phylo = TRUE,
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
  if (!is.null(tree) & phylo & !pyseer_phylo) {
    tr <- ape::keep.tip(tree, samples)
    
    if (method == "fixed") {
      D <- ape::cophenetic.phylo(tr)[samples, samples]
      dist_file <- file.path(tmpdir_time, "distances.txt")
      write.table(D, dist_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    } else if (method == "mixed") {
      D <- ape::cophenetic.phylo(tr)[samples, samples]
      sim <- exp(-D)   # similarity from distances
      sim_file <- file.path(tmpdir_time, "similarities.txt")
      write.table(sim, sim_file, sep = "\t", quote = FALSE,
                  row.names = TRUE, col.names = NA)
    } 
  } else if (!is.null(tree) & phylo & pyseer_phylo) {
    tree_file <- file.path(tmpdir_time, "core_genome.tree")
    
    # Write tree to Newick file
    ape::write.tree(tree, file = tree_file)
    
    # not usable yet...
    if (method == "fixed") {
      # compute phylo distances
      cmd_str <- paste0("python scripts/phylogeny_distance.py ", file.path(tmpdir_time, "core_genome.tree"))
      full_cmd <- paste(conda_bin, "run -n", pyseer_env, cmd_str,
                        paste0("> ", file.path(tmpdir_time, "distances.txt")))
      
      message("Running: ", full_cmd)
      system(full_cmd)
    }
    else if (method == "mixed") {
      # compute phylo distances
      cmd_str <- paste0("python scripts/phylogeny_distance.py --lmm ", file.path(tmpdir_time, "core_genome.tree"))
      full_cmd <- paste(conda_bin, "run -n", pyseer_env, cmd_str,
                        paste0("> ", file.path(tmpdir_time, "similarities.txt")))
      
      message("Running: ", full_cmd)
      system(full_cmd)
      
    }

  } else {
    if (method == "fixed") {
      # compute kinship
      K <- X %*% t(X)
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
  
  if (method == "fixed") {cmd <- c(cmd, paste0("--distances ", dist_file)) }
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
