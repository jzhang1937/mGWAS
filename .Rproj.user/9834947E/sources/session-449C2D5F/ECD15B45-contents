library(simulatr)
#source("snp.sim.mod.R")
#source("snp.sim.Q.mod.R")
#source("coalescent.sim.mod.R")
#source("methods.R")
#source("simulation_helpers.R")

# name of the simulation
sim_name <- "treeWAS_sim_settings"

# Create the parameter grid
# The default values for the simulation
seed = 1234
baseline_values = list(n.ind = 400,  # sample size
                       n.snps = 1000,
                       n.snps.assoc = 40, 
                       assoc.prob = 75)
varying_values = list(n.ind = seq(200, 600, 200),  # sample size
                      n.snps = seq(800, 1200, 200),  # dimension of Z
                      n.snps.assoc = seq(0, 60, 20),
                      assoc.prob = seq(60,90,15))
parameter_grid <- symcrt2::create_param_grid_fractional_factorial_list(
  varying_values,
  baseline_values
)

# function that inputs parameters and outputs the ground truth inferential target(s)
get_ground_truth <- function(n.snps, n.snps.assoc){
  list(nonnulls = R.utils::withSeed(sample(n.snps, n.snps.assoc), seed = seed))
}

# Add a column to parameter grid containing ground truth inferential targets
parameter_grid <- parameter_grid |> simulatr::add_ground_truth(get_ground_truth)

# The methods to compare, banded_precision gcm vs oat gcm
fixed_parameters <- list(
  B = 100,                      # number of data realizations
  seed = 4                    # seed to set prior to generating data and running methods
)

# define data-generating model 
generate_data_treeWAS_mod <- function(n.ind, n.snps, n.snps.assoc, assoc.prob, ground_truth) {
  data <- coalescent.sim.mod(n.ind = n.ind, n.snps = n.snps,
                             n.snps.assoc = n.snps.assoc,
                             assoc.prob = assoc.prob,
                             ground.truth = ground_truth)
  return(data)
}
generate_data_f <- function(n.ind, n.snps, n.snps.assoc, assoc.prob, ground_truth){
  simulate.data <- suppressWarnings(treeWAS::coalescent.sim(n.ind = n.ind, n.snps = n.snps,
                                     n.snps.assoc = n.snps.assoc,
                                     assoc.prob = assoc.prob))
  data <- list(snps = simulate.data$snps, phen = simulate.data$phen)
  data
}

# need to call simulatr_function() to give simulatr a few more pieces of info
generate_data_function <- simulatr_function(
  f = generate_data_f,                        
  arg_names = formalArgs(generate_data_f),    
  loop = TRUE
)

test.names <- c("terminal", "simultaneous", "subsequent", "cor", "fisher")
combined_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = c("terminal", "simultaneous", "subsequent", "cor", "fisher"),
                              plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

subsequent_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = "subsequent", plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

terminal_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = "terminal", plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

simultaneous_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = "simultaneous", plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

cor_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = "cor", plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

fisher_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  results <- treeWAS::treeWAS(snps = snps, phen = phen,
                              test = "fisher", plot.tree = FALSE,
                              plot.null.dist.pairs = FALSE)
  results
}

fisher_exact_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  n <- nrow(snps)
  p <- ncol(snps)
  fisher_one <- function(x, y) {
    test <- fisher.test(x = x, y = y)
    test$p.value
  }
  p.vals <- apply(X = snps, MARGIN = 2, FUN = fisher_one, y = phen)
  p.bh <- p.adjust(p = p.vals, method = "BH")
  p.bonferroni <- p.adjust(p = p.vals, method = "bonferroni")
  fisher <- list()
  fisher$p.vals <- p.vals
  fisher$bh <- which(p.bh < alpha)
  fisher$bonferroni <- which(p.bonferroni < alpha)
  results <- fisher
  results
}

chisq_f <- function(data){
  snps <- data$snps
  phen <- data$phen
  n <- nrow(snps)
  p <- ncol(snps)
  chisq_one <- function(x, y) {
    test <- chisq.test(x = x, y = y)
    test$p.value
  }
  p.vals <- apply(X = snps, MARGIN = 2, FUN = chisq_one, y = phen)
  p.bh <- p.adjust(p = p.vals, method = "BH")
  p.bonferroni <- p.adjust(p = p.vals, method = "bonferroni")
  chisq <- list()
  chisq$p.vals <- p.vals
  chisq$bh <- which(p.bh < alpha)
  chisq$bonferroni <- which(p.bonferroni < alpha)
  results <- chisq
  results
}

# create simulatr functions
combined_spec_f <- simulatr_function(f = combined_f, arg_names = character(0), loop = TRUE)
subsequent_spec_f  <- simulatr_function(f = subsequent_f, arg_names = character(0), loop = TRUE)
terminal_spec_f  <- simulatr_function(f = terminal_f, arg_names = character(0), loop = TRUE)
simultaneous_spec_f  <- simulatr_function(f = simultaneous_f, arg_names = character(0), loop = TRUE)
cor_spec_f  <- simulatr_function(f = cor_f, arg_names = character(0), loop = TRUE)
fisher_spec_f  <- simulatr_function(f = fisher_f, arg_names = character(0), loop = TRUE)
fisher_exact_spec_f  <- simulatr_function(f = fisher_exact_f, arg_names = character(0), loop = TRUE)
chisq_spec_f  <- simulatr_function(f = chisq_f, arg_names = character(0), loop = TRUE)


run_method_functions <- list(subsequent = subsequent_spec_f, 
                             terminal = terminal_spec_f,
                             simultaneous = simultaneous_spec_f,
                             cor = cor_spec_f, fisher = fisher_spec_f,
                             fisher_exact = fisher_exact_spec_f,
                             chisq = chisq_spec_f)
run_method_functions <- list(combined = combined_spec_f, 
                             fisher_exact = fisher_exact_spec_f,
                             chisq = chisq_spec_f)



evaluation_functions <- list()

simulatr_spec <- simulatr_specifier(
  parameter_grid,
  fixed_parameters,
  generate_data_function, 
  run_method_functions,
  evaluation_functions
)

check_results <- check_simulatr_specifier_object(simulatr_spec, B_in = 1)

sim_results <- check_simulatr_specifier_object(simulatr_spec)
