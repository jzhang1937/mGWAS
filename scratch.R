library(devtools)
library(treeWAS)
source("snp.sim.mod.R")
source("snp.sim.Q.mod.R")
source("coalescent.sim.mod.R")
try <- treeWAS::coalescent.sim
iter <- 100
n.ind <- 100
n.snps <- 1000
n.snps.assoc <- 400
assoc.prob <- 60
ground_truth <- list(nonnulls = R.utils::withSeed(sample(n.snps, n.snps.assoc), seed = 1))
test.names <- c("terminal", "simultaneous", "subsequent", "cor", "fisher")

for (i in 1:iter) {
  # generate tree
  tree <- treeWAS::coalescent.sim(n.ind = n.ind, n.snps = n.snps,
                                  n.snps.assoc = n.snps.assoc,
                                  assoc.prob = assoc.prob)
  # perform tests
  tests <- treeWAS(snps = tree$snps, phen = tree$phen,
                   test = test.names, plot.tree = FALSE,
                   plot.null.dist.pairs = FALSE,
                   p.value.correct = "fdr",
                   p.value.by = "count"
                   )
  
}

try <- generate_data_treeWAS_mod(n.ind = n.ind, n.snps = n.snps,
                          n.snps.assoc = n.snps.assoc,
                          assoc.prob = assoc.prob,
                          ground_truth = ground_truth)


