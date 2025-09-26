subsequent_function <- function(n.sim = 10000, p.value.correct = "fdr",
                                p.value.by = "count") {
  placeholder <- function(snps, phen) {
    test <- treeWAS::treeWAS(snps = snps, phen = phen,
                             test = "subsequent", plot.tree = FALSE,
                             plot.null.dist.pairs = FALSE,
                             p.value.correct = p.value.correct,
                             p.value.by = p.value.by)
  }
  return(placeholder)
}

terminal_function <- function(n.sim = 10000, p.value.correct = "fdr",
                              p.value.by = "count") {
  placeholder <- function(snps, phen) {
    test <- treeWAS::treeWAS(snps = snps, phen = phen,
                             test = "terminal", plot.tree = FALSE,
                             plot.null.dist.pairs = FALSE,
                             p.value.correct = p.value.correct,
                             p.value.by = p.value.by)
    test$nonnulls <- test$nonnulls$sig.snps$SNP.locus
  }
  return(placeholder)
}

simultaneous_function <- function(n.sim = 10000, p.value.correct = "fdr",
                                  p.value.by = "count") {
  placeholder <- function(snps, phen) {
    test <- treeWAS::treeWAS(snps = snps, phen = phen,
                             test = "simultaneous", plot.tree = FALSE,
                             plot.null.dist.pairs = FALSE,
                             p.value.correct = p.value.correct,
                             p.value.by = p.value.by)
  }
  return(placeholder)
}

cor_function <- function(p.value.correct = "fdr",
                         p.value.by = "count") {
  placeholder <- function(snps, phen) {
    test <- treeWAS::treeWAS(snps = snps, phen = phen,
                             test = "cor", plot.tree = FALSE,
                             plot.null.dist.pairs = FALSE,
                             p.value.correct = p.value.correct,
                             p.value.by = p.value.by)
  }
  return(placeholder)
  
}

fisher_function <- function(p.value.correct = "fdr",
                            p.value.by = "count") {
  placeholder <- function(snps, phen) {
    test <- treeWAS::treeWAS(snps = snps, phen = phen,
                             test = "fisher", plot.tree = FALSE,
                             plot.null.dist.pairs = FALSE,
                             p.value.correct = p.value.correct,
                             p.value.by = p.value.by)
  }
  return(placeholder)
}

fisher_exact <- function(snps, phen, alpha = 0.1) {
  n <- nrow(snps)
  p <- ncol(snps)
  p.vals <- apply(X = snps, MARGIN = 2, FUN = fisher_one, y = phen)
  p.bh <- p.adjust(p = p.vals, method = "BH")
  p.bonferroni <- p.adjust(p = p.vals, method = "bonferroni")
  fisher <- list()
  fisher$p.vals <- p.vals
  fisher$bh <- which(p.bh < alpha)
  fisher$bonferroni <- which(p.bonferroni < alpha)
  fisher
}

fisher_one <- function(x, y) {
  test <- fisher.test(x = x, y = y)
  test$p.value
}


chisq <- function(snps, phen, alpha = 0.1) {
  n <- nrow(snps)
  p <- ncol(snps)
  p.vals <- apply(X = snps, MARGIN = 2, FUN = chisq_one, y = phen)
  p.bh <- p.adjust(p = p.vals, method = "BH")
  p.bonferroni <- p.adjust(p = p.vals, method = "bonferroni")
  chisq <- list()
  chisq$p.vals <- p.vals
  chisq$bh <- which(p.bh < alpha)
  chisq$bonferroni <- which(p.bonferroni < alpha)
  chisq
}

chisq_one <- function(x, y) {
  test <- chisq.test(x = x, y = y)
  test$p.value
}
  
