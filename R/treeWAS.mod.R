####################
## coalescent.sim ##
####################

## a function for simulating trees under a fully-linked coalescent model.
## optional simulation of a phenotype and phenotypically-associated SNPs is implemented.
## optional use of a distribution to guide the substitution rate of the non-associated SNPs is implemented.

########################################################################

###################
## DOCUMENTATION ##
###################

#' Simulate a tree, phenotype, and genetic data.
#'
#' This funtion allows the user to simulate a phylogenetic tree, as well as
#' phenotypic and genetic data, including associated and unassociated loci.
#'
#' @param n.ind An integer specifying the number of individual genomes to simulate
#'                (ie. the number of terminal nodes in the tree).
#' @param n.snps An integer specifying the number of genetic loci to simulate.
#' @param n.subs Either an integer or a vector (containing a distribution) that is
#'                used to determine the number of substitutions
#'                to occur on the phylogenetic tree for each genetic locus (see details).
#' @param n.snps.assoc An optional integer specifying the number of genetic loci
#' @param assoc.prob An optional integer (> 0, <= 100) specifying the strength of the
#'                    association between the n.snps.assoc loci and the phenotype (see details).
#' @param n.phen.subs An integer specifying the expected number of phenotypic
#'                      substitutions to occur on the phylogenetic tree (through the same process as
#'                      the n.subs parameter when n.subs is an integer (see details)).
#' @param phen An optional vector containing a phenotype for each of the
#'              n.ind individuals if no phenotypic simulation is desired.
#' @param plot A logical indicating whether to generate a plot of the phylogenetic tree (\code{TRUE}) or not (\code{FALSE}, the default).
#' @param heatmap A logical indicating whether to produce a heatmap of the genetic distance
#'                  between the simulated genomes of the n.ind individuals.
#' @param reconstruct Either a logical indicating whether to attempt to reconstruct
#'                      a phylogenetic tree using the simulated genetic data, or one of c("UPGMA", "nj", "ml")
#'                      to specify that tree reconstruction is desired by one of these three methods
#'                      (Unweighted Pair Group Method with Arithmetic Mean, Neighbour-Joining, Maximum-Likelihood).
#' @param dist.dna.model A character string specifying the type of model to use in reconstructing the phylogenetic tree for
#'                          calculating the genetic distance between individual genomes, only used if \code{tree} is
#'                          a character string (see ?dist.dna).
#' @param grp.min An optional number between 0.1 and 0.9 to control the proportional size of the smaller phenotypic group.
#' @param row.names An optional vector containing row names for the individuals to be simulated.
#' @param set An integer (1, 2, or 3) required to select the method of generating associated loci if \code{n.snps.assoc} is not zero.
#' @param coaltree A logical indicating whether to generate a coalescent tree (\code{TRUE}, the default),
#'                 or an rtree-type tree (\code{FALSE}, see ?rtree).
#' @param s If \code{set} is 3, the \code{s} parameter controls a baseline number of substitutions to be
#'          experienced by the phenotype and associated loci: by default, 20.
#' @param af If \code{set} is 3, the \code{af} parameter provides an association factor,
#'              controlling the preference for association over non-association at associated loci:  by default, 10 (for a 10x preference).
#' @param filename.plot An optional character string denoting the file location for saving any plots produced; else \code{NULL}.
#' @param seed An optional integer to control the pseudo-randomisation process and allow for identical repeat runs of the function;
#'             else \code{NULL}.
#'
#' @details
#' \strong{Homoplasy Distribution}
#'
#' The homoplasy distribution contains the number of substitutions per site.
#'
#' If the value of the \code{n.subs} parameter is set to an integer, this integer is
#' used as the parameter of a Poisson distribution from which the number of substitutions to
#' occur on the phylogenetic tree is drawn for each of the \code{n.snps} simulated genetic loci.
#'
#' The \code{n.subs} argument can also be used to provide a distribution
#' to define the number of substitutions per site.
#'
#' It must be in the form of a \emph{named} vector (or table), or a vector in which the \emph{i}'th element
#' contains the number of \emph{loci} that have been estimated to undergo \emph{i} substitutions on the tree.
#' The vector must be of length \emph{max n.subs}, and "empty" indices must contain zeros.
#' For example: the vector \code{n.subs = c(1833, 642, 17, 6, 1, 0, 0, 1)},
#' could be used to define the homoplasy distribution for a dataset with 2500 loci,
#' where the maximum number of substitutions to be undergone on the tree by any locus is 8,
#' and no loci undergo either 6 or 7 substitutions.
#'
#'
#' \strong{Association Probability}
#'
#' The \code{assoc.prob} parameter is only functional when \code{set} is set to 1.
#' If so, \code{assoc.prob} controls the strength of association through a process analagous to dilution.
#' All \code{n.snps.assoc} loci are initially simulated to undergo a substitution
#' every time the phenotype undergoes a substitution (ie. perfect association).
#' The assoc.prob parameter then acts like a dilution factor, removing \code{(100 - assoc.prob)\%}
#' of the substitutions that occurred during simulation under perfect association.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#' @examples
#' \dontrun{
#' ## load example homoplasy distribution
#' data(dist_0)
#' str(dist_0)
#'
#' ## simulate a matrix with 10 associated loci:
#' dat <- coalescent.sim(n.ind = 100,
#'                         n.snps = 1000,
#'                         n.subs = dist_0,
#'                         n.snps.assoc = 10,
#'                         assoc.prob = 90,
#'                         n.phen.subs = 15,
#'                         phen = NULL,
#'                         plot = TRUE,
#'                         heatmap = FALSE,
#'                         reconstruct = FALSE,
#'                         dist.dna.model = "JC69",
#'                         grp.min = 0.25,
#'                         row.names = NULL,
#'                         coaltree = TRUE,
#'                         s = NULL,
#'                         af = NULL,
#'                         filename = NULL,
#'                         set = 1,
#'                         seed = 1)
#'
#' ## examine output:
#' str(dat)
#'
#' ## isolate elements of output:
#' snps <- dat$snps
#' phen <- dat$phen
#' snps.assoc <- dat$snps.assoc
#' tree <- dat$tree
#' }
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
#  @useDynLib phangorn, .registration = TRUE

############
## NOTES: ##
############
## theta_p changed to n.phen.subs (and just n.subs in phen.sim.R)


## OLD ARGS: ##
# (n.ind=100, gen.size=10000, sim.by="locus",
#  theta=1*2, dist=NULL,
#  n.phen.subs=15, phen=NULL,
#  n.snps.assoc=5, assoc.option="all", assoc.prob=90,
#  haploid=TRUE, biallelic=TRUE, seed=NULL,
#  plot=TRUE, heatmap=FALSE, plot2="UPGMA")

## NEW ARGS: ##
# n.ind <- 100
# n.snps <- 10000
# n.subs <- dist_0.1
# n.snps.assoc <- 10
# assoc.prob <- 100
# n.phen.subs <- 15
# phen <- NULL
# plot <- TRUE
# heatmap <- FALSE
# reconstruct <- FALSE
# dist.dna.model <- "JC69"
# grp.min <- 0.25
# row.names <- NULL
# coaltree <- FALSE
# s <- 20
# af <- 10
# set <- 3 # NULL #
# filename <- NULL
# seed <- 77


## TOY EG for PRESENTATION:
# c.sim <- coalescent.sim(n.ind=15,
#                         n.snps=1000,
#                         n.subs=dist_0,
#                         n.snps.assoc=10,
#                         assoc.prob=90,
#                         n.phen.subs=5,
#                         phen=NULL,
#                         plot=TRUE,
#                         heatmap=FALSE,
#                         reconstruct=FALSE,
#                         dist.dna.model="JC69",
#                         grp.min = 0.25,
#                         row.names=NULL,
#                         coaltree = TRUE,
#                         s = 10,
#                         af = 5,
#                         filename = NULL,
#                         set=1,
#                         seed=5)

## EG:
# c.sim <- coalescent.sim(n.ind=100,
#                         n.snps=10000,
#                         n.subs=dist_0,
#                         n.snps.assoc=10,
#                         assoc.prob=100,
#                         n.phen.subs=15,
#                         phen=NULL,
#                         plot=TRUE,
#                         heatmap=FALSE,
#                         reconstruct=FALSE,
#                         dist.dna.model="JC69",
#                         grp.min = 0.25,
#                         row.names=NULL,
#                         coaltree = TRUE,
#                         s = 10,
#                         af = 5,
#                         filename = NULL,
#                         set=1,
#                         seed=78)


coalescent.sim.mod <- function(n.ind = 100,
                               n.snps = 10000,
                               n.subs = 1,
                               n.snps.assoc = 0,
                               assoc.prob = 100,
                               n.phen.subs = 15,
                               phen = NULL,
                               plot = TRUE,
                               heatmap = FALSE,
                               reconstruct = FALSE,
                               dist.dna.model = "JC69",
                               grp.min = 0.25,
                               row.names = TRUE,
                               set = 1,
                               tree = NULL,
                               coaltree = TRUE,
                               s = 20,
                               af = 10,
                               filename.plot = NULL,
                               ground.truth = NULL,
                               seed = NULL){
  ## load packages:
  # require(adegenet)
  # require(ape)
  
  ## store input args to return call at end of fn:
  args <- mget(names(formals()), sys.frame(sys.nframe()))
  
  filename <- filename.plot
  
  if(is.null(set)){
    set <- 1
    # warning("set cannot be NULL. Setting set to 1.")
  }
  sets <- NULL
  
  if(length(which(c(plot, heatmap, reconstruct)==TRUE))==1){
    par(ask=FALSE)
  }else{
    par(ask=TRUE)
  }
  
  ## Allow assoc.prob to be in percent or
  ## as a proportion (-> eg 80 or 90, ie out of 100):
  if(!is.null(assoc.prob)){
    if(assoc.prob[1] >= 0 & assoc.prob <= 1){
      assoc.prob <- assoc.prob*100
    }
  }
  
  ################################
  ## Simulate Phylogenetic Tree ##
  ################################
  ## tree provided?
  if(!is.null(tree)){
    ## check:
    if(is.null(n.ind)){
      n.ind <- length(tree$tip.label)
    }
    if(length(tree$tip.label) != n.ind){
      warning("n.ind did not match length(tree$tip.label). Simulating tree instead.")
      tree <- NULL
    }
    if(!is.null(phen)){
      if(length(tree$tip.label) != length(phen)){
        warning("length(phen) did not match length(tree$tip.label). Simulating tree instead.")
        tree <- NULL
      }
    }
  }
  
  ## simulate tree:
  if(is.null(tree)){
    if(coaltree == TRUE){
      if(!is.null(seed)) set.seed(seed)
      tree <- treeWAS::coalescent.tree.sim(n.ind = n.ind, seed = seed)
    }else{
      if(!is.null(seed)) set.seed(seed)
      tree <- rtree(n = n.ind)
    }
  }
  ## Always work with tree in pruningwise order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  if(!is.rooted(tree)) tree <- midpoint(tree)
  
  ########################
  ## Simulate Phenotype ##
  ########################
  if(set == 3){
    
    ############
    ## NEW Q: ##
    ############
    ## if Q contains RATES --> P contains probs
    ## QUESTION -- HOW TO INTERPRET/PREDICT THE RELATIVE EFFECTS OF ASSOC.FACTOR AND N.SUBS (+ BRANCH LENGTH)
    ## ON ASSOC STRENGTH, N.SUBS PER TREE ?
    
    if(is.null(s)) s <- 15 # n.subs
    if(is.null(af)) af <- 10 # association factor*
    ## (*af = the relative prob of var1 changing state twd being in-phase w var2 vs. changing to the opposite state, out-of-step w var2 )
    ## Modify s to account for sum(tree$edge.length):
    ## --> ~ s/10 (= 1.5) for coaltree
    ## OR --> ~ s/100 (= 0.15) for rtree
    s <- s/sum(tree$edge.length)
    
    
    Q.mat <- matrix(c(NA,     1*s, 1*s, 0,
                      1*af*s, NA,  0,   1*af*s,
                      1*af*s, 0,   NA,  1*af*s,
                      0,      1*s, 1*s, NA),
                    nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
    
    diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
    # Q <- Q.mat
    
    ## SAVE PANEL PLOT:
    if(!is.null(filename)){
      pdf(file=filename[[2]], width=7, height=11)
      # dev.copy(pdf, file=filename[[2]], width=7, height=11)
    }
    
    ## RUN SNP.SIM.Q: ##
    if(!is.null(seed)) set.seed(seed)
    system.time(
      snps.list <- snp.sim.Q.mod(n.snps = n.snps,
                                 n.subs = n.subs,
                                 snp.root = NULL,
                                 n.snps.assoc = n.snps.assoc,
                                 assoc.prob = assoc.prob,
                                 ## dependent/corr' transition rate/prob mat:
                                 Q = Q.mat,
                                 tree = tree,
                                 n.phen.subs = n.phen.subs,
                                 phen.loci = NULL,
                                 heatmap = FALSE,
                                 reconstruct = FALSE,
                                 dist.dna.model = "JC69",
                                 grp.min = grp.min,
                                 row.names = NULL,
                                 set=set,
                                 ground.truth=ground.truth,
                                 seed=seed)
    )
    snps <- snps.list$snps
    snps.assoc <- snps.list$snps.assoc
    sets <- NULL
    phen <- snps.list$phen
    phen.nodes <- snps.list$phen.nodes
    n.mts <- snps.list$n.mts
    
    if(!is.null(filename)){
      ## end saving panel plot:
      dev.off()
    }
    
  }else{
    
    
    ##################
    ## SETS 1 and 2 ##
    ##################
    
    if(is.null(phen)){
      if(!is.null(seed)) set.seed(seed)
      ## get list of phenotype simulation output
      phen.list <- treeWAS::phen.sim(tree, n.subs = n.phen.subs, grp.min = grp.min, seed = seed)
      
      ## get phenotype for terminal nodes only
      phen <- phen.list$phen
      
      ## get phenotype for all nodes,
      ## terminal and internal
      phen.nodes <- phen.list$phen.nodes
      
      ## get the indices of phen.subs (ie. branches)
      phen.loci <- phen.list$phen.loci
    }else{
      #############################
      ## User-provided Phenotype ##
      #############################
      phen.nodes <- treeWAS::asr(phen, tree, type="parsimony")
      
      ## get COLOR for NODES
      nodeCol <- "grey"
      if(Hmisc::all.is.numeric(phen.nodes[!is.na(phen.nodes)])){
        var <- as.numeric(as.character(phen.nodes))
      }else{
        var <- as.character(phen.nodes)
      }
      levs <- unique(var[!is.na(var)])
      if(length(levs) == 2){
        ## binary:
        # myCol <- c("red", "blue")
        myCol <- c("blue", "red")
        nodeCol <- var
        ## for loop
        for(i in 1:length(levs)){
          nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
        } # end for loop
      }else{
        if(is.numeric(var)){
          ## numeric:
          myCol <- num2col(var, col.pal = seasun)
          nodeCol <- myCol
        }else{
          ## categorical...
          myCol <- funky(length(levs))
          nodeCol <- var
          ## for loop
          for(i in 1:length(levs)){
            nodeCol <- replace(nodeCol, which(nodeCol == levs[i]), myCol[i])
          } # end for loop
        }
      }
      nodeCol <- as.vector(unlist(nodeCol))
      
      ## get COLOR for EDGES
      edgeCol <- rep("black", nrow(tree$edge))
      for(i in 1:nrow(tree$edge)){
        edgeCol[i] <- nodeCol[tree$edge[i,2]]
        if(is.na(nodeCol[tree$edge[i,1]]) | is.na(nodeCol[tree$edge[i,2]])){
          edgeCol[i] <- "grey"
        }else{
          ## No grey if truly continuous...
          if(length(levs) < length(tree$tip.label)/10){
            if(nodeCol[tree$edge[i,1]] != nodeCol[tree$edge[i,2]]) edgeCol[i] <- "grey"
          }
        }
      }
      phen.loci <- which(edgeCol == "grey")
    }
    
    
    ###################
    ## Simulate SNPs ##
    ###################
    
    ## TO DO: #######################################
    ## CHECK SNP SIMULATION FOR COMPUTATIONAL SPEED
    ## 10 --> 53 --> 12.5
    ## Are the remaining extra 2.5 seconds still just a result of the while loop??
    ## Or has anything slowed down in the post-processing steps as well?
    
    #   n.snps <- 10000 # 13.3
    #   n.snps <- 100000 # 153.8
    #   n.snps <- 1000000 # >> 1941.7
    #
    
    if(!is.null(seed)) set.seed(seed)
    # system.time(
    snps.list <- snp.sim.mod(n.snps=n.snps,
                             n.subs=n.subs,
                             n.snps.assoc=n.snps.assoc,
                             assoc.prob=assoc.prob,
                             tree=tree,
                             phen.loci=phen.loci,
                             heatmap=heatmap,
                             reconstruct=reconstruct,
                             dist.dna.model=dist.dna.model,
                             row.names = NULL,
                             set=set,
                             ground.truth=ground.truth,
                             seed=seed)
    # )
    
    snps <- snps.list$snps
    snps.assoc <- snps.list$snps.assoc
    sets <- snps.list$sets
    n.mts <- snps.list$n.mts
    
  }
  
  #################################
  ## Plot Tree showing Phenotype ##
  #################################
  
  ## SAVE TREE PLOT:
  if(!is.null(filename)){
    pdf(file=filename[[1]], width=7, height=11)
    # dev.copy(pdf, file=filename[[1]], width=7, height=11)
  }
  
  if(plot==TRUE){
    if("try-error" %in% class(try(treeWAS::plot_phen(tree = tree,
                                                     phen.nodes = phen.nodes,
                                                     plot = plot)))){
      warning("Oops-- something went wrong when trying to plot
              phenotypic changes on tree.")
    }else{
      phen.plot.col <- treeWAS::plot_phen(tree = tree,
                                          phen.nodes = phen.nodes,
                                          plot = plot, main.title=FALSE)
    }
    
    ##################
    ## SET 2 CLADES ##
    ##################
    ## plot set2 clade sets along tips:
    if(!is.null(sets)){
      ## Get CLADES:
      set1 <- names(sets)[which(sets == 1)]
      set2 <- names(sets)[which(sets == 2)]
      
      tip.labs <- tree$tip.label
      if(coaltree == FALSE) tip.labs <- removeFirstN(tip.labs, 1) ## assuming tip.labs are prefaced w/ "t" for all rtrees...
      cladeCol <- rep(NA, length(tip.labs))
      cladeCol <- replace(cladeCol, which(tip.labs %in% set1), "black")
      cladeCol <- replace(cladeCol, which(tip.labs %in% set2), "grey")
      ## PLOT CLADES along tips:
      ## coaltree:
      if(coaltree == TRUE){
        tiplabels(text=NULL, cex=0.6, adj=c(0.65, 0.5), col=cladeCol, pch=15) # adj=c(0.65, 0.75)
        ## NOT SURE WHY/WHEN ADJ WORKS/w WHAT VALUES?????
      }else{
        ## rtree:
        tiplabels(text=NULL, cex=0.75, adj=c(0.65, 0.75), col=cladeCol, pch=15) # adj=c(0.65, 0.75)
        ## NOT SURE WHY/WHEN ADJ WORKS/w WHAT VALUES?????
      }
    }
    
    if(!is.null(filename)){
      ## end saving tree plot:
      dev.off()
    }
    
  }
  
  ################
  ## Handle phen
  phen.rec.method <- levs <- NULL
  phen.ori <- phen
  
  ## Convert to numeric (required for assoc tests):
  na.before <- length(which(is.na(phen)))
  
  ## CHECK for discrete vs. continuous (Only binary & continuous implemented)
  levs <- unique(as.vector(unlist(phen)))
  n.levs <- length(levs[!is.na(levs)])
  ## BINARY: ##
  if(n.levs == 2){
    ## Convert phen to numeric:
    if(!is.numeric(phen)){
      if(Hmisc::all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        phen <- as.numeric(as.factor(phen))
      }
    }
    ## Set phen.rec.method: ##
    phen.rec.method <- "discrete"
    # if(!is.null(phen.type)){
    #   if(phen.type == "continuous"){
    #     phen.rec.method <- "continuous"
    #     warning("phen is binary. Are you sure phen.type is 'continuous'?")
    #   }
    # } # end phen.type (binary)
  }else{
    ## DISCRETE or CONTINUOUS: ##
    ## Convert phen to numeric:
    if(!is.numeric(phen)){
      if(Hmisc::all.is.numeric(phen)){
        phen <- as.numeric(as.character(phen))
      }else{
        stop("phen has more than 2 levels but is not numeric (and therefore neither binary nor continuous).")
      }
    }
    ## Set phen.rec.method: ##
    phen.rec.method <- "continuous"
    ## Get proportion unique:
    prop.u <- length(unique(phen))/length(phen)
    # if(!is.null(phen.type)){
    #   if(phen.type == "discrete"){
    #     phen.rec.method <- "discrete"
    #     if(prop.u > 0.5){
    #       cat("Performing *discrete* reconstruction, although phen is ", round(prop.u, 2)*100, "% unique:
    #            Are you sure phen.type is 'discrete'?\n", sep="")
    #     }
    #   }
    # } # end phen.type (discrete/continuous)
  }
  ## ensure ind names not lost
  names(phen) <- names(phen.ori)
  
  
  ## (+) Copy treeWAS lines for all checks/cleaning: phen.rec ~1201:1244...
  ## (+) Add args: phen.reconstruction, phen.rec.method, snps.reconstruction, snps.rec.method.
  ## (+?) Add check/warning for categorical phen (if phen.type=discrete but prop.u < 0.1 (?))
  snps.rec <- treeWAS::asr(var = snps, tree = tree, unique.cols = TRUE) # , type = snps.reconstruction)
  phen.rec <- treeWAS::asr(var = phen, tree = tree, method = phen.rec.method) #, type = phen.reconstruction)
  
  ################
  
  ################
  ## Get Output ##
  ################
  out <- list(snps, snps.rec, snps.assoc, phen, phen.rec, tree, sets, args, n.mts)
  names(out) <- c("snps", "snps.rec", "snps.assoc", "phen", "phen.rec", "tree", "sets", "args", "n.mts")
  return(out)
  
} # end coalescent.sim

snp.sim.mod <- function(n.snps = 10000,
                        n.subs = 1,
                        snp.root = NULL,
                        n.snps.assoc = 0,
                        assoc.prob = 100,
                        tree = coalescent.tree.sim(100),
                        phen.loci = NULL,
                        heatmap = FALSE,
                        reconstruct = FALSE,
                        dist.dna.model = "JC69",
                        row.names = NULL,
                        set = NULL,
                        ground.truth = NULL,
                        seed = 1){
  
  # require(adegenet)
  # require(ape)
  
  ## Allow assoc.prob to be in percent or
  ## as a proportion (-> eg 80 or 90, ie out of 100):
  if(!is.null(assoc.prob)){
    if(assoc.prob[1] >= 0 & assoc.prob <= 1){
      assoc.prob <- assoc.prob*100
    }
  }
  
  ##################
  ## HANDLE TREE: ##
  ##################
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)
  
  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs (if binary tree), from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
  ####################################################################
  
  
  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  n.ind <- min(tree$edge[,1])-1 # tree$Nnode+1
  gen.size <- n.snps
  edges <- tree$edge
  
  if(!is.null(seed)) set.seed(seed)
  
  ## Simulate genotype for root individual: ##
  
  ## For n.subs = n or = dist approaches:
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  
  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)
  
  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA
  
  snps.assoc <- NULL
  
  ## if n.snps.assoc is neither NULL nor 0:
  if(is.null(n.snps.assoc)) n.snps.assoc <- 0
  if(n.snps.assoc != 0){
    
    ## get non.assoc gen.size
    gen.size.ori <- gen.size
    gen.size <- gen.size-n.snps.assoc
    
    ## assign snps.assoc to be the last n.snps.assoc snps columns
    snps.assoc <- c((gen.size+1):(gen.size+n.snps.assoc))
  }
  
  ###################
  ## Handle n.subs ##
  ###################
  
  ## Either an integer
  ## --> draw n.subs from a Poisson distribution w parameter n.subs
  ## OR a vector (containing a distribution)
  ## --> use this distribution to define n.subs-per-site
  
  if(length(n.subs)==1 & is.null(names(n.subs))){
    
    #####################
    ## NO DISTRIBUTION ##
    #####################
    ## if no distribution is inputted,
    ## use normal simulation procedure
    ## (ie. Poisson parameter 1):
    
    warning("Using n.subs as Poisson parameter because input n.subs was of length 1 and had no names.")
    
    ## draw the number of mutations to occur at each site:
    if(!is.null(seed)) set.seed(seed)
    n.mts <- rpois(n=gen.size, lambda=(n.subs))
    ## for any n.mts==0, re-sample
    if(any(n.mts == 0)){
      ## need to change seed or we'll get trapped in the while loop
      if(!is.null(seed)) seed.i <- seed
      for(i in 1:length(n.mts)){
        while(n.mts[i]==0){
          if(!is.null(seed)){
            seed.i <- seed.i+1
            set.seed(seed.i)
          }
          n.mts[i] <- rpois(n=1, lambda=(n.subs))
        }
      }
    }
    
  }else{
    
    ###############################################
    ## DISTRIBUTION or RATES (fitPagel Q matrix) ##
    ###############################################
    
    ##################
    ## DISTRIBUTION ##
    ##################
    
    ## if a distribution is provided by the user,
    ## we use this to determine the number of substitutions
    ## to occur at what proportion of sites (note that
    ## we may not be simulating the same number of sites)
    
    dist <- n.subs
    
    ## check for names first!
    if(!is.null(names(dist))){
      ## only modify if names are numeric
      if(all.is.numeric(names(dist))){
        noms <- as.numeric(names(dist))
        aligned <- sapply(c(1:length(dist)), function(e) noms[e] == e)
        ## if any names do not correspond to their index, add zeros where missing:
        if(any(aligned == FALSE)){
          dist.new <- rep(0, max(noms))
          dist.new[noms] <- dist
          names(dist.new) <- c(1:length(dist.new))
          dist <- dist.new
        }
      }
    } # end check for missing places
    
    ## get dist.prop, a distribution containing the counts
    ## of the number of SNPs to be simulated that will have
    ## i many substitutions
    dist.sum <- sum(dist)
    dist.prop <- round((dist/dist.sum)*gen.size)
    ## check that these counts sum to gen.size,
    ## else add the remainder to the largest n.subs count
    ## (OR should we just add these to the n.subs=1 set ??? ###
    ## likely to be the same thing, but could not be...)
    if(sum(dist.prop) != gen.size){
      m <- which.max(dist.prop)
      #m <- 1
      if(sum(dist.prop) < gen.size){
        dist.prop[m] <- dist.prop[m] + (gen.size - sum(dist.prop))
      }
      if(sum(dist.prop) > gen.size){
        dist.prop[m] <- dist.prop[m] - (sum(dist.prop) - gen.size)
      }
    }
    
    ## get rid of useless trailing 0s
    while(dist.prop[length(dist.prop)] == 0){
      dist.prop <- dist.prop[c(1:(length(dist.prop)-1))]
    }
    
    ## make n.mts, a vector of length ncol(snps)
    n.mts <- rep(1111, gen.size)
    loci.available <- c(1:gen.size)
    ## assign dist.prop[i] elements of n.mts
    ## to be the same as the n.subs
    ## indicated by i, the given element of dist.prop
    for(j in 1:length(dist.prop)){
      ## provided there are not 0 sites to have this number of substitutions...
      if(dist.prop[j] > 0){
        if(length(loci.available) > 1){
          if(!is.null(seed)) set.seed(seed)
          ## assign dist.prop[i] elements of n.mts to be i
          loci.selected <- sample(loci.available, dist.prop[j], replace = FALSE)
          loci.available <- loci.available[-which(loci.available %in% loci.selected)]
        }else{
          ## if there is only 1 (the last) loci available,
          ## we select this one:
          loci.selected <- loci.available
        }
        n.mts[loci.selected] <- j
      }
    }
    ## Remove unnecessary objects...
    rm(dist)
    rm(dist.prop)
    rm(dist.sum)
    # } # end dist
  } # end fitPagel or dist
  
  ############################
  ## Assign mts to branches ##
  ############################
  
  if(n.snps.assoc != 0){
    if(length(phen.loci) > 0){
      ## for snps.assoc (the last n.snps.assoc snps, for now),
      ## add n.mts == n.phen.loci s.t these sites mutate at each
      ## and every phen.loci (for now, to be diluted later
      ## according to assoc.prob if !=100)
      n.mts <- c(n.mts, rep(length(phen.loci), n.snps.assoc))
    }else{
      stop("No phen.loci specified, so no associated loci can be generated.")
    }
  }
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #############################################################################
  ## GENERATE ALL SNPs FIRST, THEN REPLACE ANY NON-POLYMORPHIC IN WHILE LOOP ##
  #############################################################################
  
  #############################
  ## GET NON-ASSOCIATED SNPS ##
  #############################
  
  ## for each site, draw the branches to which
  ## you will assign the mts for this site
  ## (~ branch length):
  
  l.edge <- length(tree$edge.length)
  ## Get vector of FALSEs of length tree$edge.length:
  null.vect <- rep(FALSE, l.edge)
  
  
  ## TO DO: Memory inefficient step... Improve if possible?
  ## NB: this condition should never be met within treeWAS | parsimonious homoplasy dist:
  if(max(n.mts) > l.edge){
    if(!is.null(seed)) set.seed(seed)
    snps.loci <- NULL
    for(e in 1:length(n.mts)){
      snps.loci <- list()
      repTF <- FALSE
      if(n.mts[e] > l.edge) repTF <- TRUE
      snps.loci[[e]] <- replace(null.vect, sample(c(1:l.edge),
                                                  n.mts[e],
                                                  replace=repTF,
                                                  prob=tree$edge.length), TRUE)
      
      snps.loci <- t(do.call(rbind, snps.loci))
    }
    
  }else{
    if(!is.null(seed)) set.seed(seed)
    snps.loci <- sapply(c(1:length(n.mts)),
                        function(e)
                          replace(null.vect,
                                  sample(c(1:l.edge),
                                         n.mts[e],
                                         replace=FALSE,
                                         prob=tree$edge.length), TRUE))
  }
  
  ## rearrange snps.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## locations of the mutations that will
  ## occur on that branch
  snps.loci <- sapply(c(1:nrow(snps.loci)),
                      function(e) which(snps.loci[e,] == TRUE))
  
  
  ## get the node names for all individuals (terminal and internal)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  # we will store the output in a list called snps:
  snps <- list()
  ## we start w all inds having same genotype as root:
  snps[all.inds] <- rep(list(gen.root), length(all.inds))
  
  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  snps.loci.ori <- snps.loci
  ## will need to treat repeat loci differently...
  snps.loci.unique <- lapply(snps.loci, unique)
  
  
  #############################
  ## For Loop to get new nts ##
  #############################
  for(i in x){
    ## for all snps other than root, we mutate the
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!treeWAS::.is.integer0(snps.loci.unique[[i]])){
      new.nts[[i]] <- !snps[[tree$edge[i,1]]][snps.loci.unique[[i]]]
      
      
      ## if any loci are selected for multiple mutations
      ## within their given branch length:
      if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
        ## identify which loci are repeaters
        repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
        ## how many times they repeat
        n.reps <- repeats - 1
        ## the positions of these loci in the vector of snps loci
        toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
        ## run chain of re-sampling to end in our new nt for repeater loci:
        foo <- list()
        for(j in 1:length(toRepeat)){
          foo[[j]] <- new.nts[[i]][toRepeat[j]]
          for(k in 1:n.reps[j]){
            if(k==1){
              foo[[j]][k] <- !foo[[j]][1]
              
            }else{
              foo[[j]][k] <- !foo[[j]][k-1]
            }
          }
          ## retain only the last nt selected
          out <- sapply(c(1:length(foo)),
                        function(e) foo[[e]][length(foo[[e]])])
        }
        ## for the loci with repeated mts, replace these positions
        ## in new.nts with the corresponding elements of out, above.
        new.nts[[i]][toRepeat] <- out
      } # end of if statement for repeaters
      
      ## update ancestral genotype with new.nts:
      temp <- snps[[tree$edge[i,1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]
      snps[[tree$edge[i,2]]] <- temp
      
    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      snps[[tree$edge[i,2]]] <- snps[[tree$edge[i,1]]]
    }
  } # end of for loop selecting new nts at mutator loci
  
  ####################################################
  ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres) ##
  ####################################################
  
  ## temporarily assemble non-associated loci into matrix:
  # temp.ori <- do.call("rbind", snps)
  temp <- do.call("rbind", snps)
  
  ## keep only rows containing terminal individuals:
  # temp.ori <- temp.ori[1:n.ind, ]
  temp <- temp[1:n.ind, ]
  
  ######################################################################################################################################################################
  
  ## identify n.minor.allele required to meet polyThres:
  polyThres <- 0.01
  n.min <- n.ind*polyThres
  
  ## make a list of any NON-polymorphic loci:
  csum <- colSums(temp)
  toRepeat <- which(csum < n.min | csum > (nrow(temp) - n.min))
  
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #################################################################
  ## REPLACE ANY NON-POLYMORPHIC LOCI & GENERATE ASSOCIATED SNPS ##
  #################################################################
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  
  #################################################
  ## REPLACE NON-POLYMORPHIC NON-ASSOCIATED SNPS ##
  #################################################
  
  ## just working w rows containing inds (not internal nodes)
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  ####################
  ## NEW while loop ##
  ####################
  
  counter <- 0
  if(!is.null(seed)) seed.i <- seed
  
  while(length(toRepeat) > 0){
    # print("LENGTH toREPEAT"); print(length(toRepeat))
    ## for each site, draw the branches to which
    ## you will assign the mts for this site
    ## (~ branch length):
    
    ## Get vector of FALSEs of length tree$edge.length:
    null.vect <- rep(FALSE, length(tree$edge.length))
    
    
    if(max(n.mts[toRepeat]) > length(tree$edge.length)){
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=TRUE,
                                           prob=tree$edge.length), TRUE))
    }else{
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=FALSE,
                                           prob=tree$edge.length), TRUE))
    }
    
    ## rearrange snps.loci s.t it becomes a
    ## list of length tree$edge.length,
    ## each element of which contains the
    ## locations of the mutations that will
    ## occur on that branch
    snps.loci <- sapply(c(1:nrow(snps.loci)),
                        function(e) which(snps.loci[e,] == TRUE))
    
    
    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))
    # we will store the output in a list called snps:
    # snps <- list()
    ## we start w all inds having same genotype as root:
    # snps[all.inds][toRepeat] <- rep(list(gen.root[toRepeat]), length(all.inds))
    for(i in 1:length(all.inds)){
      snps[[all.inds[i]]][toRepeat] <- gen.root[toRepeat]
    }
    
    ## store replacement nts in list new.nts:
    new.nts <- list()
    ## distinguish btw list of loci and unique list
    snps.loci.ori <- snps.loci
    ## will need to treat repeat loci differently...
    snps.loci.unique <- lapply(snps.loci, unique) # (identical to snps.loci?)
    
    
    #############################
    ## For Loop to get new nts ##
    #############################
    for(i in x){
      ## for all snps other than root, we mutate the
      ## genome of the node preceding it, according to snps.loci.
      ## Draw new nts for each locus selected for mutation:
      if(!treeWAS::.is.integer0(snps.loci.unique[[i]])){
        new.nts[[i]] <- !snps[[tree$edge[i,1]]][toRepeat][snps.loci.unique[[i]]]
        
        ## if any loci are selected for multiple mutations
        ## within their given branch length:
        if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
          ## identify which loci are repeaters
          repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
          ## how many times they repeat
          n.reps <- repeats - 1
          ## the positions of these loci in the vector of snps loci
          toRep <- which(snps.loci.unique[[i]] %in% names(repeats))
          ## run chain of re-sampling to end in our new nt for repeater loci:
          foo <- list()
          for(j in 1:length(toRep)){
            foo[[j]] <- new.nts[[i]][toRep[j]]
            for(k in 1:n.reps[j]){
              if(k==1){
                foo[[j]][k] <- !foo[[j]][1]
              }else{
                foo[[j]][k] <- !foo[[j]][k-1]
              }
            }
            ## retain only the last nt selected
            out <- sapply(c(1:length(foo)),
                          function(e) foo[[e]][length(foo[[e]])])
          }
          ## for the loci with repeated mts, replace these positions
          ## in new.nts with the corresponding elements of out, above.
          new.nts[[i]][toRep] <- out
        } # end of if statement for repeaters
        
        ## update ancestral genotype with new.nts:
        toto <- snps[[tree$edge[i,1]]][toRepeat]
        toto[snps.loci.unique[[i]]] <- new.nts[[i]]
        snps[[tree$edge[i,2]]][toRepeat] <- toto
        
      }else{
        ## if no mts occur on branch, set genotype of
        ## downstream individual to be equal to ancestor's
        snps[[tree$edge[i,2]]][toRepeat] <- snps[[tree$edge[i,1]]][toRepeat]
      }
    } # end of for loop selecting new nts at mutator loci
    
    
    ## temporarily assemble non-associated loci into matrix:
    # temp.new <- do.call("rbind", snps)
    temp <- do.call("rbind", snps)
    
    ## keep only rows containing terminal individuals:
    # temp <- temp.new[1:n.ind, ]
    temp <- temp[1:n.ind, ]
    
    ##############################################################
    # temp <- temp.new
    
    ######################################
    ##### while loop CHECK here: #########
    ######################################
    ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres)
    
    ## identify n.minor.allele required to meet polyThres:
    polyThres <- 0.01
    n.min <- n.ind*polyThres
    
    ## make a list of any NON-polymorphic loci:
    toRepeat.ori <- toRepeat
    temp.toRepeat <- temp[, toRepeat.ori]
    
    ## make a vector of any NON-polymorphic loci:
    ## If ncol = 1:
    if(!is.matrix(temp.toRepeat)){
      csum <- sum(temp.toRepeat)
      if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
        toRepeat <- toRepeat.ori
      }else{
        toRepeat <- NULL
      }
    }else{
      ## if temp.toRepeat is a true matrix:
      if(ncol(temp.toRepeat) > 0){
        csum <- colSums(temp.toRepeat)
        toRepeat <- toRepeat.ori[which(csum < n.min | csum > (nrow(temp.toRepeat) - n.min))]
      }else{
        csum <- sum(temp.toRepeat)
        if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
          toRepeat <- toRepeat.ori
        }else{
          toRepeat <- NULL
        }
      }
    }
    
    counter <- counter+1
    # print("COUNTER"); print(counter)
    # print("toRepeat"); print(length(toRepeat))
    
    
  } # end NEW while loop...
  #######################################
  ### while loop ENDS here: #############
  #######################################
  
  colnames(temp) <- c(1:ncol(temp))
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  
  #########################
  ## GET ASSOCIATED SNPS ##
  #########################
  
  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)
  if(n.snps.assoc != 0){
    ## get snps.loci for the ASSOCIATED snps (ie. set to phen.loci) ##
    for(i in 1:n.snps.assoc){
      ## recall: phen.loci contains the tree EDGES on which phen subs occur
      subs.edges <- phen.loci
      
      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]
      
      ## get nt for each individual at this locus
      ## assign to (and replace) the snps.assoc elements of loci
      temp[, snps.assoc[i]] <- treeWAS:::.get.locus(subs.edges = subs.edges,
                                                    root.nt = root.nt,
                                                    tree = tree)[1:n.ind]
    }
  } # end of snps.assoc generation
  
  
  ###########################################
  ## GET COMPLETE SNPS MATRIX ("snps"): ##
  ###########################################
  
  ## Remove unnecessary objects...
  rm(snps)
  ## Create snps matrix:
  snps <- temp
  ## Remove unnecessary objects...
  rm(temp)
  
  
  ## keep only rows containing terminal individuals:
  snps <- snps[1:n.ind, ]
  
  ## clear memory
  # gc()
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  
  
  ###############################################
  ## MODIFY SNPS.ASSOC ACCORDING TO ASSOC.PROB ##
  ###############################################
  
  if(n.snps.assoc != 0){
    ## if we have any imperfect associations... ##
    if(any(assoc.prob != 100)){
      ## check length
      if(length(assoc.prob) != n.snps.assoc){
        ## if only 1 prob value given...
        if(length(assoc.prob) == 1){
          ## ... assume uniform assoc.prob;
          # assoc.prob <- rep(assoc.prob, n.snps.assoc)
          ## ... OR, randomly draw assoc.probs from normal dist around assoc.prob;
          assoc.prob <- round(rnorm(n.snps.assoc, assoc.prob, 2), 3)
          ## no warning needed
        }else{
          ## BUT if assoc.prob of random length:
          ## repeat until of length n.snps.assoc
          assoc.prob <- rep(assoc.prob, length.out=n.snps.assoc)
          ## and print warning (only if not of length n.snps.assoc OR 1)
          warning("assoc.prob not of length n.snps.assoc;
                  sequence will be repeated until correct length is reached.")
        }
      } # end checks
      
      
      ## re-pseudo-randomise seed:
      if(!is.null(seed)){
        seed.i <- seed*c(1:n.snps.assoc)*10
      }
      
      ## for each associated SNP,
      ## we undo some associations | assoc.prob for that snp.assoc
      for(i in 1:n.snps.assoc){
        
        ## re-pseudo-randomise seed:
        if(!is.null(seed)){
          set.seed(seed.i[i])
        }
        
        prob <- assoc.prob[i]
        ## only if the association is imperfect
        if(prob <= 100){
          
          # print("prob"); print(prob)
          # print("inv. prob"); print((1 - (prob/100)))
          
          ## draw snps to change at snps.assoc[i]
          n.toChange <- round(nrow(snps)*(1 - (prob/100)))
          
          # print("N to change"); print(n.toChange)
          # print("Nrow SNPS"); print(str(snps))
          
          toChange <- sample(c(1:nrow(snps)), n.toChange)  # , replace=TRUE
          
          ## change those snps at rows toChange, loci snps.assoc[i]
          for(j in 1:length(toChange)){
            snps[toChange[j], snps.assoc[i]] <-
              treeWAS::selectBiallelicSNP(snps[toChange[j], snps.assoc[i]])
          } # end for loop
        }
      } # end for loop
    } # end any assoc.prob != 100
  } # end modification | assoc.prob
  
  
  ##############################
  ## PLOTS & TREECONSTRUCTION ##
  ##############################
  tree.reconstructed <- NULL
  if(heatmap == TRUE || reconstruct!=FALSE){
    
    ## CONVERT TO CHARACTER: ##
    dna <- snps
    dna <- replace(dna, which(dna == TRUE), "a")
    dna <- replace(dna, which(dna == "FALSE"), "t")
    
    ## Get DNAbin object:
    dna <- as.DNAbin(dna)
    # rownames(dna) <- c(1:nrow(snps))
    
    
    #############
    ## HEATMAP ##
    #############
    if(heatmap==TRUE){
      heatmap.DNAbin(dna=dna,
                     dist.dna.model=dist.dna.model)
    }
    
    ##########################################
    ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
    ##########################################
    if(reconstruct!=FALSE){
      tree.reconstructed <- tree.reconstruct(dna, # [1:n.ind,]
                                             method=reconstruct,
                                             dist.dna.model=dist.dna.model,
                                             plot=TRUE)
    }
    
    ## Remove unnecessary object:
    rm(dna)
    
  } # end heatmap, treeconstruction
  
  ##################
  ## CONVERT SNPS ##
  ##################
  
  ## (NO LONGER) CONVERT TO NUMERIC: ##
  ## Convert from logical to binary SNPs (for terminal nodes only):
  # snps <- replace(snps, which(snps == TRUE), 1)
  
  ## Reassort snps.assoc to new columns:
  if(!is.null(snps.assoc)){
    
    ## update snps.assoc to reflect true loci
    gen.size.final <- ncol(snps)
    snps.assoc.loci.ori <- c((gen.size.final-(n.snps.assoc-1)):gen.size.final)
    
    #########################################
    ## RANDOMIZE SNPS.ASSOC LOCI POSITIONS ##
    #########################################
    
    ## Re-enabled snps.assoc loci "randomization" by
    ## just drawing indices and shuffling the columns accordingly...
    ## draw which SNPs will be associated to the phenotype
    if(!is.null(seed)) set.seed(seed)
    
    if(is.null(ground.truth)) {
      snps.assoc.loci <- sort(sample(c(1:gen.size.final),
                                     n.snps.assoc,
                                     replace=FALSE))
    } else {
      snps.assoc.loci <- sort(ground.truth$nonnulls)
    }
    
    
    
    snps.indices <- c(1:gen.size.final)
    snps.ori <- snps
    
    snps.non.assoc <- snps[,c(1:(gen.size.final-n.snps.assoc))]
    snps.assoc <- snps[,snps.assoc.loci.ori]
    snps.new <- matrix(99, nrow=nrow(snps), ncol=gen.size.final)
    snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
    snps.new[,snps.assoc.loci] <- snps.assoc
    snps <- snps.new
    snps.assoc <- snps.assoc.loci
    
  } # end snps.assoc randomization
  
  ###############################
  ## Assign row & column names ##
  ###############################
  
  ## assign/generate row.names
  if(!is.null(row.names)){
    ## If row.names have been provided in args list, assign them:
    if(length(row.names) == nrow(snps)){
      ## match tree$tip.label?
      if(!is.null(tree$tip.label)){
        if(all(row.names %in% tree$tip.label) & all(tree$tip.label %in% row.names)){
          ## REORDER to match tree$tip.labs if possible:
          if(!identical(row.names, tree$tip.label)){
            ord <- match(tree$tip.label, row.names)
            row.names <- row.names[ord]
          }
        }
      }
      rownames(snps) <- row.names
    }else{
      if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
    }
  }else{
    ## Otherwise, try to assign rownames(snps) to match tree$tip.label:
    if(!is.null(tree$tip.label)){
      if(length(tree$tip.label) == nrow(snps)){
        rownames(snps) <- tree$tip.label
      }else{
        rownames(snps) <- 1:nrow(snps)
        warning("The length of tree$tip.label was not equal to nrow(snps) being simulated;
              rownames(snps) have been set to 1:N and will not match tree$tip.label.")
      }
    }
  }
  
  ## generate column names:
  colnames(snps) <- 1:ncol(snps)
  
  
  #################################################
  ## SIM SET 2 (complementary clade-wise assoc): ##
  #################################################
  sets <- NULL
  if(!is.null(snps.assoc)){
    if(!is.null(set)){
      if(set == 2){
        
        ## Want to divide tree into 2 sets of clades btw 1/3:2/3 and 1/2:1/2
        clades <- tab <- grp.options <- sets.complete <- list()
        
        min.size <- ceiling((n.ind)*(1/3))
        max.size <- floor((n.ind)*(2/3))
        grp1 <- n.ind
        
        ## get 2 sets of clades:
        
        ##########################################
        ## Pick SETS for ANY TREE (pretty sure) ##
        ##########################################
        dec <- grp <- sets.temp <- sets.complete <- list()
        
        inds <- c(1:(n.ind))
        # new.root <- tree$edge[1,1] # initial root
        new.root <- n.ind+1 # initial root
        
        counter <- 0
        #######################################
        ## WHILE LOOP to get size of clades: ##
        #######################################
        while(grp1 < min.size | grp1 > max.size){
          
          ## get all descendants of root node:
          all.dec <- treeWAS:::.getDescendants(tree, node=new.root)
          
          ## get all descendants in first 2 major clades:
          dec[[1]] <- treeWAS:::.getDescendants(tree, node=all.dec[1])
          dec[[2]] <- treeWAS:::.getDescendants(tree, node=all.dec[2])
          
          ## get terminal inds only:
          sets.temp[[1]] <- dec[[1]][which(dec[[1]] %in% inds)]
          sets.temp[[2]] <- dec[[2]][which(dec[[2]] %in% inds)]
          
          grp[[1]] <- length(sets.temp[[1]])
          grp[[2]] <- length(sets.temp[[2]])
          
          max.grp <- which.max(c(grp[[1]], grp[[2]]))
          new.root <- all.dec[max.grp]
          
          set1 <- sets.temp[[max.grp]]
          
          sets <- rep(2, length(inds))
          sets <- replace(sets, set1, 1)
          names(sets) <- tree$tip.label
          
          counter <- counter+1
          
          grp1 <- grp[[max.grp]]
          
        } # end while loop
        
        
        set1 <- names(sets)[which(sets == 1)]
        set2 <- names(sets)[which(sets == 2)]
        ###########
        
        ########################
        ## MODIFY SNPS.ASSOC: ##
        ########################
        snps.assoc.set1 <- 1:round(length(snps.assoc)/2)
        snps.assoc.set2 <- (round(length(snps.assoc)/2)+1):length(snps.assoc)
        
        ## replace set1 snps with 0 at all inds in clade.set1:
        for(e in 1:length(snps.assoc.set1)){
          # snps[which(rownames(snps) %in% set1), snps.assoc[snps.assoc.set1[e]]] <- 0
          snps[which(rownames(snps) %in% set1), snps.assoc[snps.assoc.set1[e]]] <- FALSE
        }
        ## replace set2 snps with 0 at all inds in clade.set2:
        for(e in 1:length(snps.assoc.set2)){
          # snps[which(rownames(snps) %in% set2), snps.assoc[snps.assoc.set2[e]]] <- 0
          snps[which(rownames(snps) %in% set2), snps.assoc[snps.assoc.set2[e]]] <- FALSE
        }
      }
    }
  } # end sim set 2
  
  
  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, tree.reconstructed, sets, n.mts)
  names(out) <- c("snps", "snps.assoc", "tree.reconstructed", "sets", "n.mts")
  
  return(out)
  
} # end snp.sim

## WARNING:
## SOME LINES OF SNP.SIM.Q.R HAVE FALLEN OUT OF DATE w SNP.SIM.R.
## BEFORE USING SNP.SIM.Q, UPDATE (at least!):
## sapply, add n.mts > l.edge check,
## convert snps to logical, for loop (selectBiallelicSNP --> !l[[i]])

###############
## snp.sim.Q ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Aternative SNPs simulation fn.
#'
#' Currently under development. Please use the regular snp.sim function to simulate genetic data.
#'
#' @param n.snps An integer specifying the number of snps columns to be simulated.
#'               At present, all snps are simulated to be correlated to the
#'               provided \code{phen.reconstruction} along the provided \code{tree}.
#' @param tree A \code{phylo} object containing the phylogenetic tree; or, a character string,
#'                one of \code{"NJ"}, \code{"BIONJ"} (the default), or \code{"parsimony"};
#'                or, if NAs are present in the distance matrix, one of: \code{"NJ*"} or \code{"BIONJ*"},
#'                specifying the method of phylogenetic reconstruction.
#' @param phen.reconstruction Either a character string specifying \code{"parsimony"} (the default) or \code{"ML"} (maximum likelihood)
#'                              for the ancestral state reconstruction of the phenotypic variable,
#'                              or a vector containing this reconstruction if it has been performed elsewhere.
#' @param s If \code{set} is 3, the \code{s} parameter controls a baseline number of substitutions to be
#'          experienced by the phenotype and associated loci: by default, 20.
#' @param af If \code{set} is 3, the \code{af} parameter provides an association factor,
#'              controlling the preference for association over non-association at associated loci:  by default, 10 (for a 10x preference).
#' @param plot A logical indicating whether to generate a plot of the phylogenetic tree (\code{TRUE}) or not (\code{FALSE}, the default).
#' @param heatmap A logical indicating whether to produce a heatmap of the genetic distance
#'                  between the simulated genomes of the n.ind individuals.
#' @param reconstruct Either a logical indicating whether to attempt to reconstruct
#'                      a phylogenetic tree using the simulated genetic data, or one of c("UPGMA", "nj", "ml")
#'                      to specify that tree reconstruction is desired by one of these three methods
#'                      (Unweighted Pair Group Method with Arithmetic Mean, Neighbour-Joining, Maximum-Likelihood).
#' @param dist.dna.model A character string specifying the type of model to use in reconstructing the phylogenetic tree for
#'                          calculating the genetic distance between individual genomes, only used if \code{tree} is
#'                          a character string (see ?dist.dna).
#' @param grp.min (Not yet (re-)implemented in this function.)
#'                An optional number between 0.1 and 0.9 to control the proportional size of the smaller phenotypic group.
#' @param row.names An optional vector containing row names for the individuals to be simulated.
#' @param seed An optional integer to control the pseudo-randomisation process and allow for identical repeat runs of the function;
#'             else \code{NULL}.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @examples
#' ## Example ##
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
#  @useDynLib phangorn, .registration = TRUE

snp.sim.Q.mod.new <- function(n.snps = 10,
                              tree = coalescent.tree.sim(100, seed=1),
                              phen.reconstruction, ## provide phen.rec and use this to set s=n.phen.subs...
                              s = NULL, ## n.subs for correlated snps.assoc (s over-rides n.phen.subs)
                              af = 10, ## association factor (relative odds of assoc vs. non-assoc in Q mat)
                              snp.root = NULL,
                              heatmap = FALSE,
                              reconstruct = FALSE,
                              dist.dna.model = "JC69",
                              grp.min = 0.25,
                              row.names = NULL,
                              ground.truth = NULL,
                              seed=1){
  
  # require(adegenet)
  # require(ape)
  
  temp <- snps <- snps.assoc <- tree.reconstructed <- sets <- phen <- phen.nodes <- NULL
  
  n.snps.assoc <- n.snps
  
  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)
  
  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
  ####################################################################
  
  
  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  n.ind <- min(tree$edge[,1])-1 # tree$Nnode+1
  gen.size <- n.snps
  edges <- tree$edge
  
  if(!is.null(seed)) set.seed(seed)
  
  ## Simulate genotype (& phen) for root individual: ##
  
  ## if snp.root given:
  if(!is.null(snp.root)){
    if(length(snp.root) == 1){
      ## select only root state --> different SNP sim method?
      if(snp.root %in% c(0, FALSE)) gen.root <- rep(FALSE, gen.size)
      if(snp.root %in% c(1, TRUE)) gen.root <- rep(TRUE, gen.size)
    }else{
      ## if snp.root provided for all loci:
      if(length(snp.root) == gen.size){
        if(length(unique(snp.root[!is.na(snp.root)])) == 2){
          gen.root <- as.logical(as.numeric(as.factor(snp.root))-1)
        }else{
          warning("snp.root must be binary; ignoring.")
        }
      }else{
        warning("snp.root should either be of length 1 or length n.snps; ignoring.")
      }
    }
  }
  
  ## check phen.rec:
  if(!is.null(phen.reconstruction)){
    if(length(unique(phen.reconstruction[!is.na(phen.reconstruction)])) == 2){
      p.noms <- names(phen.reconstruction)
      phen.reconstruction <- as.logical(as.numeric(phen.reconstruction))
      names(phen.reconstruction) <- p.noms
      
      ## Infer n.phen.subs:
      n.phen.subs <- length(which(phen.reconstruction[tree$edge[,1]] != phen.reconstruction[tree$edge[,2]]))
    }else{
      phen.reconstruction <- NULL
      n.phen.subs <- NULL
      warning("phen.reconstruction must be binary. Ignorning (sorry).")
    }
  }else{
    phen.reconstruction <- NULL
    n.phen.subs <- NULL
  }
  
  
  ## For n.subs = n or = dist approaches:
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  if(is.null(phen.reconstruction)){
    phen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  }else{
    phen.root <- phen.reconstruction[min(tree$edge[,1])]
  }
  
  
  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)
  
  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA
  
  snps.assoc <- NULL
  
  ## if n.snps.assoc is neither NULL nor 0:
  if(is.null(n.snps.assoc)) n.snps.assoc <- 0
  if(n.snps.assoc != 0){
    
    ## get non.assoc gen.size
    gen.size.ori <- gen.size
    gen.size <- gen.size-n.snps.assoc
    
    ## assign snps.assoc to be the last n.snps.assoc snps columns
    snps.assoc <- c((gen.size+1):(gen.size+n.snps.assoc))
  }
  
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  
  #########################
  ## SIM ASSOCIATED SNPS ##
  #########################
  
  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)
  
  snps.assoc.nodes <- phen.nodes <- NULL
  
  if(n.snps.assoc != 0){
    
    ## (*) Now assuming phen.rec provided as input
    ## If s provided, s overrides n.phen.subs;
    ## else s=n.phen.subs as inferred from phen.reconstruction:
    if(is.null(s)){
      if(!is.null(n.phen.subs)){
        s <- n.phen.subs
      }else{
        s <- 15
        warning("Setting s (n.subs for snps.assoc in Q rate matrix) as 15.
                User must provide phen.reconstruction (s=n.phen.subs) or set s argument to specify manually.")
      }
    }
    ## Scale s to account for total edge.length:
    s <- s/sum(tree$edge.length)
    
    ## Set af = the association factor:
    ## *af = the relative odds of var1 changing state towards being in-phase w var2 vs.
    ## changing to the opposite state, out-of-step w var2, over a given branch length...
    ## (Q is actually an instantaneous transition rate matrix, but we will account for branch lengths later.)
    if(is.null(af)){
      af <- 10
      warning("Setting association factor af = 10.")
    }
    
    
    ## Get Q, the dependent transition rate/prob matrix:
    Q.mat <- matrix(c(NA,     1*s, 1*s, 0,
                      1*af*s, NA,  0,   1*af*s,
                      1*af*s, 0,   NA,  1*af*s,
                      0,      1*s, 1*s, NA),
                    nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
    
    diag(Q.mat) <- sapply(c(1:nrow(Q.mat)), function(e) -sum(Q.mat[e, c(1:ncol(Q.mat))[-e]]))
    Q <- Q.mat
    
    
    
    ## get snps.loci for all ASSOCIATED snps, conditional on phen, Q: ##
    snps.assoc.nodes <- list()
    N.OVERLAP <- list()
    
    Qt <- list()
    
    ## get snps.loci for FIRST ASSSOCIATED snp WITH PHEN.loci:
    i <- 1
    edges <- tree$edge
    root.nt <- gen.root[snps.assoc[i]]
    
    
    ## Need to update tree$edge.length probs to mirror sampling without replacement...
    ## BUT only remove branches from the calculation once they have had a sub on them...?
    # N.SUBS.TOTAL <- rpois(1, n.phen.subs)  ############
    
    ###############
    ## FOR LOOP: ##
    ###############
    set.seed(seed)
    for(i in 1:n.snps.assoc){
      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]
      
      snp.node <- as.list(rep(NA, length(unique(as.vector(edges)))))
      names(snp.node) <- names(phen.reconstruction)
      
      snp.node[[edges[x[1], 1]]] <- root.nt
      
      ## go from last to first edge in edges:
      for(e in x){
        
        probs <- NULL
        
        ####################################################
        ## get conditional probs for each edge w matexpo! ##
        ####################################################
        Qt[[e]] <- matexpo(Q*tree$edge.length[e])
        
        P <- Qt[[e]]
        rownames(P) <- rownames(Qt[[e]]) <- rownames(Q)
        colnames(P) <- colnames(Qt[[e]]) <- colnames(Q)
        
        if(snp.node[[edges[e, 1]]] == FALSE & phen.reconstruction[[edges[e, 1]]] == FALSE) probs <- P[1,]
        if(snp.node[[edges[e, 1]]] == FALSE & phen.reconstruction[[edges[e, 1]]] == TRUE) probs <- P[2,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.reconstruction[[edges[e, 1]]] == FALSE) probs <- P[3,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.reconstruction[[edges[e, 1]]] == TRUE) probs <- P[4,]
        
        ## Now we KNOW, a priori, the phen.node for the DESCENDANT
        probs.mod <- replace(probs, which(!as.logical(as.numeric(treeWAS::keepLastN(colnames(Q), 1))) == phen.reconstruction[[edges[e, 2]]]), 0)
        SP.dec <- sample(colnames(Q), 1, prob = probs.mod)
        
        S.dec <- as.logical(as.numeric(treeWAS::keepFirstN(SP.dec, 1)))
        names(S.dec) <- NULL
        
        snp.node[[edges[e, 2]]] <- S.dec
        
      } # end for (e) loop
      
      ## STORE SNPS.ASSOC (FOR ALL NODES):
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))
      
      
      ################################################
      ##### ##### ##### ##### ##### ##### ##### #####
      ## Get proportion overlap btw phen and snps.assoc.i:
      N.overlap.pos <- length(which(phen.reconstruction[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      # N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.overlap.neg <- n.ind - N.overlap.pos
      if(N.overlap.pos > N.overlap.neg){
        N.overlap <- +round((N.overlap.pos/n.ind), 2)
      }else{
        N.overlap <- -round((N.overlap.neg/n.ind), 2)
      }
      
      N.OVERLAP[[i]] <- N.overlap
      
      ##### ##### ##### ##### ##### ##### ##### #####
      ################################################
      
    } # end for loop
    
    ## Bind SNPs.ASSOC into matrix:
    snps.assoc.nodes <- do.call("cbind", snps.assoc.nodes)
    rownames(snps.assoc.nodes) <- names(phen.reconstruction)
    
    
    ###################################################
    ## TEMP -- COMPARE PHEN & ALL SNPS.ASSOC w PLOT: ##
    ###################################################
    par(mfrow=c(2,6))
    plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
    if(n.snps.assoc >= 10){
      ## just plot first 10:
      for(i in 1:5){
        plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE, show.axis=FALSE,
                  main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=2.8, font.main=1, cex.main=0.8)
      }
      plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
      for(i in 6:10){
        plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE, show.axis=FALSE,
                  main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=2.8, font.main=1, cex.main=0.8)
      }
    }else{
      ## plot up to 5:
      if(n.snps.assoc <= 5){
        par(mfrow=c(1,(n.snps.assoc+1)))
        for(i in 1:n.snps.assoc){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE, show.axis=FALSE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=2.8, font.main=1)
        }
      }else{
        ## plot btw 5 and 10:
        for(i in 1:5){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE, show.axis=FALSE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=2.8, font.main=1)
        }
        plot_phen(tree, phen.nodes=phen.reconstruction, main.title="phen")
        for(i in 6:n.snps.assoc){
          plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE, show.axis=FALSE,
                    main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=2.8, font.main=1)
        }
      }
    }
    
    
    par(mfrow=c(1,1)) # end temp panel plot
    
    gc()
    
  } # end of snps.assoc generation
  
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  
  
  ###########################################
  ## GET COMPLETE SNPS MATRIX ("genomes"): ##
  ###########################################
  
  if(!is.null(temp)){
    ## Create snps matrix:
    snps <- temp
    rm(temp)
    
    ## Attach snps.assoc loci to last column:
    if(!is.null(snps.assoc.nodes)){
      snps <- cbind(snps[,c(1:(ncol(snps)-(n.snps.assoc)))], snps.assoc.nodes[1:n.ind,])
    }
    ## keep only rows containing terminal individuals:
    snps <- snps[1:n.ind, ]
  }else{
    
    ## If only snps.assoc simulated (no non-assoc snps):
    snps <- snps.assoc.nodes[1:n.ind,]
  }
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  
  
  ##############################
  ## PLOTS & TREECONSTRUCTION ##
  ##############################
  tree.reconstructed <- NULL
  if(heatmap == TRUE || reconstruct!=FALSE){
    
    ## CONVERT TO CHARACTER: ##
    dna <- snps
    dna <- replace(dna, which(dna == TRUE), "a")
    dna <- replace(dna, which(dna == "FALSE"), "t")
    
    ## Get DNAbin object:
    dna <- as.DNAbin(dna)
    # rownames(dna) <- c(1:nrow(snps))
    
    
    #############
    ## HEATMAP ##
    #############
    if(heatmap==TRUE){
      heatmap.DNAbin(dna=dna,
                     dist.dna.model=dist.dna.model)
    }
    
    ##########################################
    ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
    ##########################################
    if(reconstruct!=FALSE){
      tree.reconstructed <- tree.reconstruct(dna, # [1:n.ind,]
                                             method=reconstruct,
                                             dist.dna.model=dist.dna.model,
                                             plot=TRUE)
    }
    
    ## Remove unnecessary object:
    rm(dna)
    
  } # end heatmap, treeconstruction
  
  
  ##################
  ## CONVERT SNPS ##
  ##################
  
  ## (NO LONGER) CONVERT TO NUMERIC: ##
  ## Convert from logical to binary SNPs (for terminal nodes only):
  # snps <- replace(snps, which(snps == TRUE), 1)
  
  ## Reassort snps.assoc to new columns:
  if(!is.null(snps.assoc)){
    
    if(!is.null(temp)){
      rm(temp)
      
      ## update snps.assoc to reflect true loci
      gen.size.final <- ncol(snps)
      snps.assoc.loci.ori <- c((gen.size.final-(n.snps.assoc-1)):gen.size.final)
      
      #########################################
      ## RANDOMIZE SNPS.ASSOC LOCI POSITIONS ##
      #########################################
      
      ## Re-enabled snps.assoc loci "randomization" by
      ## just drawing indices and shuffling the columns accordingly...
      ## draw which SNPs will be associated to the phenotype
      if(!is.null(seed)) set.seed(seed)
      if (is.null(ground.truth)) {
        snps.assoc.loci <- sort(sample(c(1:gen.size.final),
                                       n.snps.assoc,
                                       replace=FALSE))
      } else {
        snps.assoc.loci <- sort(ground.truth$nonnulls)
      }
      
      
      snps.indices <- c(1:gen.size.final)
      snps.ori <- snps
      
      snps.non.assoc <- snps[,c(1:(gen.size.final-n.snps.assoc))]
      snps.assoc <- snps[,snps.assoc.loci.ori]
      snps.new <- matrix(99, nrow=nrow(snps), ncol=gen.size.final)
      snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
      snps.new[,snps.assoc.loci] <- snps.assoc
      snps <- snps.new
      snps.assoc <- snps.assoc.loci
    }
  } # end snps.assoc randomization
  
  ###############################
  ## Assign row & column names ##
  ###############################
  
  ## assign/generate row.names
  if(!is.null(row.names)){
    ## If row.names have been provided in args list, assign them:
    if(length(row.names) == nrow(snps)){
      ## match tree$tip.label?
      if(!is.null(tree$tip.label)){
        if(all(row.names %in% tree$tip.label) & all(tree$tip.label %in% row.names)){
          ## REORDER to match tree$tip.labs if possible:
          if(!identical(row.names, tree$tip.label)){
            ord <- match(tree$tip.label, row.names)
            row.names <- row.names[ord]
          }
        }
      }
      rownames(snps) <- row.names
    }else{
      if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
    }
  }else{
    ## Otherwise, try to assign rownames(snps) to match tree$tip.label:
    if(!is.null(tree$tip.label)){
      if(length(tree$tip.label) == nrow(snps)){
        rownames(snps) <- tree$tip.label
      }else{
        rownames(snps) <- 1:nrow(snps)
        warning("The length of tree$tip.label was not equal to nrow(snps) being simulated;
              rownames(snps) have been set to 1:N and will not match tree$tip.label.")
      }
    }
  }
  
  ## generate column names:
  colnames(snps) <- 1:ncol(snps)
  
  ## Get association +/-:
  snps.assoc.direction <- NULL
  if(!is.null(N.OVERLAP)){
    snps.assoc.direction <- as.vector(unlist(N.OVERLAP))
  }
  
  
  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, snps.assoc.nodes, snps.assoc.direction, tree.reconstructed, phen, phen.nodes, n.mts)
  names(out) <- c("snps", "snps.assoc", "snps.assoc.nodes", "snps.assoc.direction", "tree.reconstructed",  "phen", "phen.nodes", "n.mts")
  
  return(out)
  
} # end snp.sim.Q

## WARNING:
## SOME LINES OF SNP.SIM.Q.R HAVE FALLEN OUT OF DATE w SNP.SIM.R.
## BEFORE USING SNP.SIM.Q, UPDATE (at least!):
## sapply, add n.mts > l.edge check,
## convert snps to logical, for loop (selectBiallelicSNP --> !l[[i]])

###############
## snp.sim.Q ##
###############

########################################################################

###################
## DOCUMENTATION ##
###################

#' Aternative SNPs simulation fn.
#'
#' NOT currently in use. Please use the regular snp.sim function to simulate genetic data.
#'
#' @param n.snps An integer specifying the number of snps columns to be simulated.
#'
#' @author Caitlin Collins \email{caitiecollins@@gmail.com}
#'
#'
#' @examples
#' ## Example ##
#'
#' @import adegenet
#' @rawNamespace import(ape, except = zoom)
#' @importFrom Hmisc all.is.numeric
#' @importFrom phangorn midpoint
#'
#' @export

########################################################################
#  @useDynLib phangorn, .registration = TRUE

## ARGUMENTS: ##

## n.subs <- either an integer or a vector containing a distribution of n.subs-per-site
## phen.loci <- a vector containing the indices of the edges on which phen subs occurred

#########
# tree <- get(load("/media/caitiecollins/88CC9BCECC9BB4C2/Cait 2016/Work/Xavier/Sims/set1/set1_31_tree.Rdata"))

# ## ARGS (Eg): ##
# n.snps = 1000
# n.subs = 1
# snp.root = NULL
# n.snps.assoc = 10
# assoc.prob = 100
# # ## dependent/corr' transition rate/prob mat:
# # # Q = matrix(c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2),
# # #            nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2))
# ### tree = coalescent.tree.sim(100)
# phen.loci = NULL
# n.phen.subs <- 15
# heatmap = FALSE
# reconstruct = FALSE
# dist.dna.model = "JC69"
# grp.min <- 0.25
# row.names = NULL
# set=3
# seed=1

snp.sim.Q.mod <- function(n.snps = 10000,
                          n.subs = 1,
                          snp.root = NULL,
                          n.snps.assoc = 10,
                          assoc.prob = 100,
                          ## dependent/correlated transition rate/prob mat:
                          Q = matrix(c(2, 0.75, 0.75, 1, 3, 0.5, 0.25,  3, 3, 0.25, 0.5, 3, 1, 0.75, 0.75, 2),
                                     nrow=4, byrow=T, dimnames=rep(list(c("0|0", "0|1", "1|0", "1|1")), 2)),
                          tree = coalescent.tree.sim(100),
                          n.phen.subs = 15,
                          phen.loci = NULL,
                          heatmap = FALSE,
                          reconstruct = FALSE,
                          dist.dna.model = "JC69",
                          grp.min = 0.25,
                          row.names = NULL,
                          set=3,
                          ground.truth = NULL,
                          seed=1){
  
  # require(adegenet)
  # require(ape)
  
  
  ## HANDLE TREE: ##
  ## Always work with trees in "pruningwise" order:
  tree <- reorder.phylo(tree, order="pruningwise")
  ## Trees must be rooted:
  # if(!is.rooted(tree)) tree <- midpoint(tree)
  
  ####################################################################
  ############################
  ## Get Anc-Des EDGE ORDER ##
  ############################
  ## Get sequence from lowest ("root", Nterm+1) to highest ancestral node:
  ix <- c(min(tree$edge[,1]):max(tree$edge[,1]))
  ## Get for loop index of rows in tree$edge[,1], in pairs, from lowest to highest:
  x <- as.vector(unlist(sapply(c(1:length(ix)), function(e) which(tree$edge[,1] == ix[e]))))
  ####################################################################
  
  
  ##################################
  ## GET MUTATIONS' branch & loci ##
  ##################################
  n.ind <- min(tree$edge[,1])-1 # tree$Nnode+1
  gen.size <- n.snps
  edges <- tree$edge
  
  if(!is.null(seed)) set.seed(seed)
  
  ## Simulate genotype (& phen) for root individual: ##
  
  ## if snp.root given:
  if(!is.null(snp.root)){
    if(length(snp.root) == 1){
      ## select only root state --> different SNP sim method (???)
      if(snp.root %in% c(0, FALSE)) gen.root <- rep(FALSE, gen.size)
      if(snp.root %in% c(1, TRUE)) gen.root <- rep(TRUE, gen.size)
    }else{
      ## if snp.root provided for all loci:
      if(length(snp.root) == gen.size){
        if(length(unique(snp.root[!is.na(snp.root)])) == 2){
          gen.root <- as.logical(as.numeric(as.factor(snp.root))-1)
        }else{
          warning("snp.root must be binary; ignoring.")
        }
      }else{
        warning("snp.root should either be of length 1 or length n.snps; ignoring.")
      }
    }
  }
  
  ## For n.subs = n or = dist approaches:
  gen.root <- sample(c(TRUE, FALSE), gen.size, replace=TRUE)
  phen.root <- sample(c(TRUE, FALSE), 1)
  
  ## get the sum of all branch lengths in the tree:
  time.total <- sum(tree$edge.length)
  
  ## make dummy variables in which to store the resulting n.mts variables:
  L <- lambda <- n.mts <- new.nt <- NA
  
  snps.assoc <- NULL
  
  ## if n.snps.assoc is neither NULL nor 0:
  if(is.null(n.snps.assoc)) n.snps.assoc <- 0
  if(n.snps.assoc != 0){
    
    ## get non.assoc gen.size
    gen.size.ori <- gen.size
    gen.size <- gen.size-n.snps.assoc
    
    ## assign snps.assoc to be the last n.snps.assoc snps columns
    snps.assoc <- c((gen.size+1):(gen.size+n.snps.assoc))
  }
  
  ###################
  ## Handle n.subs ##
  ###################
  
  ## Either an integer
  ## --> draw n.subs from a Poisson distribution w parameter n.subs
  ## OR a vector (containing a distribution)
  ## --> use this distribution to define n.subs-per-site
  
  if(length(n.subs)==1 & is.null(names(n.subs))){
    
    #####################
    ## NO DISTRIBUTION ##
    #####################
    ## if no distribution is inputted,
    ## use normal simulation procedure
    ## (ie. Poisson parameter 1):
    
    warning("Using n.subs as Poisson parameter because input n.subs was of length 1 and had no names.")
    
    ## draw the number of mutations to occur at each site:
    if(!is.null(seed)) set.seed(seed)
    n.mts <- rpois(n=gen.size, lambda=(n.subs))
    ## for any n.mts==0, re-sample
    if(any(n.mts == 0)){
      ## need to change seed or we'll get trapped in the while loop
      if(!is.null(seed)) seed.i <- seed
      for(i in 1:length(n.mts)){
        while(n.mts[i]==0){
          if(!is.null(seed)){
            seed.i <- seed.i+1
            set.seed(seed.i)
          }
          n.mts[i] <- rpois(n=1, lambda=(n.subs))
        }
      }
    }
    
  }else{
    
    ###############################################
    ## DISTRIBUTION or RATES (fitPagel Q matrix) ##
    ###############################################
    
    ##################
    ## DISTRIBUTION ##
    ##################
    
    ## if a distribution is provided by the user,
    ## we use this to determine the number of substitutions
    ## to occur at what proportion of sites (note that
    ## we may not be simulating the same number of sites)
    
    dist <- n.subs
    
    ## check for names first!
    if(!is.null(names(dist))){
      ## only modify if names are numeric
      if(all.is.numeric(names(dist))){
        noms <- as.numeric(names(dist))
        aligned <- sapply(c(1:length(dist)), function(e) noms[e] == e)
        ## if any names do not correspond to their index, add zeros where missing:
        if(any(aligned == FALSE)){
          dist.new <- rep(0, max(noms))
          dist.new[noms] <- dist
          names(dist.new) <- c(1:length(dist.new))
          dist <- dist.new
        }
      }
    } # end check for missing places
    
    ## get dist.prop, a distribution containing the counts
    ## of the number of SNPs to be simulated that will have
    ## i many substitutions
    dist.sum <- sum(dist)
    dist.prop <- round((dist/dist.sum)*gen.size)
    ## check that these counts sum to gen.size,
    ## else add the remainder to the largest n.subs count
    ## (OR should we just add these to the n.subs=1 set ??? ###
    ## likely to be the same thing, but could not be...)
    if(sum(dist.prop) != gen.size){
      m <- which.max(dist.prop)
      #m <- 1
      if(sum(dist.prop) < gen.size){
        dist.prop[m] <- dist.prop[m] + (gen.size - sum(dist.prop))
      }
      if(sum(dist.prop) > gen.size){
        dist.prop[m] <- dist.prop[m] - (sum(dist.prop) - gen.size)
      }
    }
    
    ## get rid of useless trailing 0s
    while(dist.prop[length(dist.prop)] == 0){
      dist.prop <- dist.prop[c(1:(length(dist.prop)-1))]
    }
    
    ## make n.mts, a vector of length ncol(snps)
    n.mts <- rep(1111, gen.size)
    loci.available <- c(1:gen.size)
    ## assign dist.prop[i] elements of n.mts
    ## to be the same as the n.subs
    ## indicated by i, the given element of dist.prop
    for(j in 1:length(dist.prop)){
      ## provided there are not 0 sites to have this number of substitutions...
      if(dist.prop[j] > 0){
        if(length(loci.available) > 1){
          if(!is.null(seed)) set.seed(seed)
          ## assign dist.prop[i] elements of n.mts to be i
          loci.selected <- sample(loci.available, dist.prop[j], replace = FALSE)
          loci.available <- loci.available[-which(loci.available %in% loci.selected)]
        }else{
          ## if there is only 1 (the last) loci available,
          ## we select this one:
          loci.selected <- loci.available
        }
        n.mts[loci.selected] <- j
      }
    }
    ## Remove unnecessary objects...
    rm(dist)
    rm(dist.prop)
    rm(dist.sum)
    # } # end dist
  } # end fitPagel or dist
  
  ############################
  ## Assign mts to branches ##
  ############################
  
  if(n.snps.assoc != 0){
    ## for snps.assoc (the last n.snps.assoc snps, for now),
    ## add n.mts == n.phen.loci s.t these sites mutate at each
    ## and every phen.loci (for now, to be diluted later
    ## according to assoc.prob if !=100)
    n.mts <- c(n.mts, rep(1, n.snps.assoc))
    # n.mts <- c(n.mts, rep(length(phen.loci), n.snps.assoc))
  }
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #############################################################################
  ## GENERATE ALL SNPs FIRST, THEN REPLACE ANY NON-POLYMORPHIC IN WHILE LOOP ##
  #############################################################################
  
  #############################
  ## GET NON-ASSOCIATED SNPS ##
  #############################
  
  ## for each site, draw the branches to which
  ## you will assign the mts for this site
  ## (~ branch length):
  
  l.edge <- length(tree$edge.length)
  ## Get vector of FALSEs of length tree$edge.length:
  null.vect <- rep(FALSE, l.edge)
  
  
  ## TO DO: Memory inefficient step... Improve if possible?
  if(max(n.mts) > l.edge){
    if(!is.null(seed)) set.seed(seed)
    snps.loci <- NULL
    for(e in 1:length(n.mts)){
      snps.loci <- list()
      repTF <- FALSE
      if(n.mts[e] > l.edge) repTF <- TRUE
      snps.loci[[e]] <- replace(null.vect, sample(c(1:l.edge),
                                                  n.mts[e],
                                                  replace=repTF,
                                                  prob=tree$edge.length), TRUE)
      
      snps.loci <- t(do.call(rbind, snps.loci))
    }
    
  }else{
    if(!is.null(seed)) set.seed(seed)
    snps.loci <- sapply(c(1:length(n.mts)),
                        function(e)
                          replace(null.vect,
                                  sample(c(1:l.edge),
                                         n.mts[e],
                                         replace=FALSE,
                                         prob=tree$edge.length), TRUE))
  }
  
  ## rearrange snps.loci s.t it becomes a
  ## list of length tree$edge.length,
  ## each element of which contains the
  ## locations of the mutations that will
  ## occur on that branch
  snps.loci <- sapply(c(1:nrow(snps.loci)),
                      function(e) which(snps.loci[e,] == TRUE))
  
  
  ## get the node names for all individuals (terminal and internal)
  all.inds <- sort(unique(as.vector(unlist(tree$edge))))
  # we will store the output in a list called snps:
  snps <- list()
  ## we start w all inds having same genotype as root:
  snps[all.inds] <- rep(list(gen.root), length(all.inds))
  
  ## store replacement nts in list new.nts:
  new.nts <- list()
  ## distinguish btw list of loci and unique list
  snps.loci.ori <- snps.loci
  ## will need to treat repeat loci differently...
  snps.loci.unique <- lapply(snps.loci, unique)
  
  
  #############################
  ## For Loop to get new nts ##
  #############################
  for(i in x){
    ## for all snps other than root, we mutate the
    ## genome of the node preceding it, according to snps.loci.
    ## Draw new nts for each locus selected for mutation:
    if(!treeWAS::.is.integer0(snps.loci.unique[[i]])){
      new.nts[[i]] <- !snps[[tree$edge[i,1]]][snps.loci.unique[[i]]]
      
      
      ## if any loci are selected for multiple mutations
      ## within their given branch length:
      if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
        ## identify which loci are repeaters
        repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
        ## how many times they repeat
        n.reps <- repeats - 1
        ## the positions of these loci in the vector of snps loci
        toRepeat <- which(snps.loci.unique[[i]] %in% names(repeats))
        ## run chain of re-sampling to end in our new nt for repeater loci:
        foo <- list()
        for(j in 1:length(toRepeat)){
          foo[[j]] <- new.nts[[i]][toRepeat[j]]
          for(k in 1:n.reps[j]){
            if(k==1){
              foo[[j]][k] <- !foo[[j]][1]
              
            }else{
              foo[[j]][k] <- !foo[[j]][k-1]
            }
          }
          ## retain only the last nt selected
          out <- sapply(c(1:length(foo)),
                        function(e) foo[[e]][length(foo[[e]])])
        }
        ## for the loci with repeated mts, replace these positions
        ## in new.nts with the corresponding elements of out, above.
        new.nts[[i]][toRepeat] <- out
      } # end of if statement for repeaters
      
      ## update ancestral genotype with new.nts:
      temp <- snps[[tree$edge[i,1]]]
      temp[snps.loci.unique[[i]]] <- new.nts[[i]]
      snps[[tree$edge[i,2]]] <- temp
      
    }else{
      ## if no mts occur on branch, set genotype of
      ## downstream individual to be equal to ancestor's
      snps[[tree$edge[i,2]]] <- snps[[tree$edge[i,1]]]
    }
  } # end of for loop selecting new nts at mutator loci
  
  ####################################################
  ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres) ##
  ####################################################
  
  ## temporarily assemble non-associated loci into matrix:
  # temp.ori <- do.call("rbind", snps)
  temp <- do.call("rbind", snps)
  
  ## keep only rows containing terminal individuals:
  # temp.ori <- temp.ori[1:n.ind, ]
  temp <- temp[1:n.ind, ]
  
  ######################################################################################################################################################################
  
  ## identify n.minor.allele required to meet polyThres:
  polyThres <- 0.01
  n.min <- n.ind*polyThres
  
  ## make a list of any NON-polymorphic loci:
  csum <- colSums(temp)
  toRepeat <- which(csum < n.min | csum > (nrow(temp) - n.min))
  
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  #################################################################
  ## REPLACE ANY NON-POLYMORPHIC LOCI & GENERATE ASSOCIATED SNPS ##
  #################################################################
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  
  #################################################
  ## REPLACE NON-POLYMORPHIC NON-ASSOCIATED SNPS ##
  #################################################
  
  ## just working w rows containing inds (not internal nodes)
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  ####################
  ## NEW while loop ##
  ####################
  
  counter <- 0
  if(!is.null(seed)) seed.i <- seed
  
  while(length(toRepeat) > 0){
    
    # print("LENGTH toREPEAT"); print(length(toRepeat))
    # print("toRepeat NAs?"); print(length(which(is.na(toRepeat))))
    # print("n.mts NAs?"); print(length(which(is.na(n.mts))))
    # print("str(n.mts)"); print(str(n.mts))
    # print("n.mts[toRepeat]"); print(n.mts[toRepeat])
    # print("max(n.mts[toRepeat])"); print(max(n.mts[toRepeat]))
    # print("length(tree$edge.length)"); print(length(tree$edge.length))
    
    ## for each site, draw the branches to which
    ## you will assign the mts for this site
    ## (~ branch length):
    
    ## Get vector of FALSEs of length tree$edge.length:
    null.vect <- rep(FALSE, length(tree$edge.length))
    
    
    if(max(n.mts[toRepeat]) > length(tree$edge.length)){
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=TRUE,
                                           prob=tree$edge.length), TRUE))
    }else{
      if(!is.null(seed)){
        seed.i <- seed.i+1
        set.seed(seed.i)
      }
      snps.loci <- sapply(c(1:length(n.mts[toRepeat])),
                          function(e)
                            replace(null.vect,
                                    sample(c(1:length(tree$edge.length)),
                                           n.mts[toRepeat][e],
                                           replace=FALSE,
                                           prob=tree$edge.length), TRUE))
    }
    
    ## rearrange snps.loci s.t it becomes a
    ## list of length tree$edge.length,
    ## each element of which contains the
    ## locations of the mutations that will
    ## occur on that branch
    snps.loci <- sapply(c(1:nrow(snps.loci)),
                        function(e) which(snps.loci[e,] == TRUE))
    
    
    ## get the node names for all individuals (terminal and internal)
    all.inds <- sort(unique(as.vector(unlist(tree$edge))))
    # we will store the output in a list called snps:
    # snps <- list()
    ## we start w all inds having same genotype as root:
    # snps[all.inds][toRepeat] <- rep(list(gen.root[toRepeat]), length(all.inds))
    for(i in 1:length(all.inds)){
      snps[[all.inds[i]]][toRepeat] <- gen.root[toRepeat]
    }
    
    ## store replacement nts in list new.nts:
    new.nts <- list()
    ## distinguish btw list of loci and unique list
    snps.loci.ori <- snps.loci
    ## will need to treat repeat loci differently...
    snps.loci.unique <- lapply(snps.loci, unique) # (identical to snps.loci?)
    
    
    #############################
    ## For Loop to get new nts ##
    #############################
    for(i in x){
      ## for all snps other than root, we mutate the
      ## genome of the node preceding it, according to snps.loci.
      ## Draw new nts for each locus selected for mutation:
      if(!treeWAS::.is.integer0(snps.loci.unique[[i]])){
        new.nts[[i]] <- !snps[[tree$edge[i,1]]][toRepeat][snps.loci.unique[[i]]]
        
        ## if any loci are selected for multiple mutations
        ## within their given branch length:
        if(length(snps.loci.ori[[i]]) != length(snps.loci.unique[[i]])){
          ## identify which loci are repeaters
          repeats <- table(snps.loci.ori[[i]])[which(table(snps.loci.ori[[i]])!=1)]
          ## how many times they repeat
          n.reps <- repeats - 1
          ## the positions of these loci in the vector of snps loci
          toRep <- which(snps.loci.unique[[i]] %in% names(repeats))
          ## run chain of re-sampling to end in our new nt for repeater loci:
          foo <- list()
          for(j in 1:length(toRep)){
            foo[[j]] <- new.nts[[i]][toRep[j]]
            for(k in 1:n.reps[j]){
              if(k==1){
                foo[[j]][k] <- !foo[[j]][1]
              }else{
                foo[[j]][k] <- !foo[[j]][k-1]
              }
            }
            ## retain only the last nt selected
            out <- sapply(c(1:length(foo)),
                          function(e) foo[[e]][length(foo[[e]])])
          }
          ## for the loci with repeated mts, replace these positions
          ## in new.nts with the corresponding elements of out, above.
          new.nts[[i]][toRep] <- out
        } # end of if statement for repeaters
        
        ## update ancestral genotype with new.nts:
        toto <- snps[[tree$edge[i,1]]][toRepeat]
        toto[snps.loci.unique[[i]]] <- new.nts[[i]]
        snps[[tree$edge[i,2]]][toRepeat] <- toto
        
      }else{
        ## if no mts occur on branch, set genotype of
        ## downstream individual to be equal to ancestor's
        snps[[tree$edge[i,2]]][toRepeat] <- snps[[tree$edge[i,1]]][toRepeat]
      }
    } # end of for loop selecting new nts at mutator loci
    
    
    ## temporarily assemble non-associated loci into matrix:
    # temp.new <- do.call("rbind", snps)
    temp <- do.call("rbind", snps)
    
    ## keep only rows containing terminal individuals:
    # temp <- temp.new[1:n.ind, ]
    temp <- temp[1:n.ind, ]
    
    ##############################################################
    # temp <- temp.new
    
    ######################################
    ##### while loop CHECK here: #########
    ######################################
    ## CHECK IF ALL LOCI ARE POLYMORPHIC (|polyThres)
    
    ## identify n.minor.allele required to meet polyThres:
    polyThres <- 0.01
    n.min <- n.ind*polyThres
    
    ## make a list of any NON-polymorphic loci:
    toRepeat.ori <- toRepeat
    temp.toRepeat <- temp[, toRepeat.ori]
    
    ## make a vector of any NON-polymorphic loci:
    ## If ncol = 1:
    if(!is.matrix(temp.toRepeat)){
      csum <- sum(temp.toRepeat)
      if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
        toRepeat <- toRepeat.ori
      }else{
        toRepeat <- NULL
      }
    }else{
      ## if temp.toRepeat is a true matrix:
      if(ncol(temp.toRepeat) > 0){
        csum <- colSums(temp.toRepeat)
        toRepeat <- toRepeat.ori[which(csum < n.min | csum > (nrow(temp.toRepeat) - n.min))]
      }else{
        csum <- sum(temp.toRepeat)
        if(csum < n.min | csum > (length(temp.toRepeat) - n.min)){
          toRepeat <- toRepeat.ori
        }else{
          toRepeat <- NULL
        }
      }
    }
    
    counter <- counter+1
    # print("COUNTER"); print(counter)
    # print("toRepeat"); print(length(toRepeat))
    
    
  } # end NEW while loop...
  #######################################
  ### while loop ENDS here: #############
  #######################################
  
  colnames(temp) <- c(1:ncol(temp))
  
  # print("WHILE LOOP DONE; toRepeat?"); print(length(toRepeat))
  ######################################################################################################################################################################
  ######################################################################################################################################################################
  
  
  
  #########################
  ## GET ASSOCIATED SNPS ##
  #########################
  
  ## Need to treat ASSOCIATED SNPs differently:
  ## (non.assoc.snps do NOT need to pass the "while" check;
  ## they just need to match phen.loci at this point.)
  
  snps.assoc.nodes <- phen.nodes <- NULL
  
  if(n.snps.assoc != 0){
    
    snps.assoc.nodes <- list()
    N.OVERLAP <- list()
    
    Qt <- list()
    
    ## get snps.loci for FIRST ASSSOCIATED snp WITH PHEN.loci:
    i <- 1
    
    edges <- tree$edge
    
    root.nt <- gen.root[snps.assoc[i]]
    
    
    ## Need to update tree$edge.length probs to mirror sampling without replacement...
    ## BUT only remove branches from the calculation once they have had a sub on them...?
    # N.SUBS.TOTAL <- 15
    N.SUBS.TOTAL <- rpois(1, n.phen.subs)  ############
    
    ################
    ## WHILE LOOP ##
    ################
    toRepeat <- TRUE
    ## WHILE LOOP STARTS HERE:
    while(toRepeat){ ############################
      
      snp.node <- phen.node <- as.list(rep(NA, length(unique(as.vector(edges)))))
      
      snp.node[[edges[x[1], 1]]] <- root.nt
      phen.node[[edges[x[1], 1]]] <- phen.root
      
      # probs.e <- round((tree$edge.length/sum(tree$edge.length)),3)
      # probs.e <- (tree$edge.length/sum(tree$edge.length))
      probs.e <- NULL
      
      N.SUBS.COUNTER <- 0
      
      ##################
      ## FOR (e) LOOP ##
      ##################
      ## go from last to first edge in edges:
      for(e in x){
        
        # print("E"); print(e)
        
        probs <- NULL
        toKeep <- toSub <- NULL
        
        ## WANT -- 15 subs/tree (thus 183 no-sub branches): ##
        
        ####################################################
        ## get conditional probs for each edge w matexpo! ##
        ####################################################
        ## (run within code...)
        Qt[[e]] <- matexpo(Q*tree$edge.length[e])
        
        P <- Qt[[e]]
        rownames(P) <- rownames(Qt[[e]]) <- rownames(Q)
        colnames(P) <- colnames(Qt[[e]]) <- colnames(Q)
        
        ##############################
        ## SNP.anc = 0 Phen.anc = 0 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == FALSE){
          probs <- P[1,]
        }
        ##############################
        ## SNP.anc = 0 Phen.anc = 1 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == TRUE){
          probs <- P[2,]
        }
        ##############################
        ## SNP.anc = 1 Phen.anc = 0 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == FALSE){
          probs <- P[3,]
        }
        ##############################
        ## SNP.anc = 1 Phen.anc = 1 ##
        ##############################
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == TRUE){
          probs <- P[4,]
        }
        
        SP.dec <- sample(colnames(Q), 1, prob = probs)
        
        S.dec <- as.logical(as.numeric(treeWAS::keepFirstN(SP.dec, 1)))
        P.dec <- as.logical(as.numeric(treeWAS::keepLastN(SP.dec, 1)))
        names(S.dec) <- names(P.dec) <- NULL
        
        snp.node[[edges[e, 2]]] <- S.dec
        phen.node[[edges[e, 2]]] <- P.dec
        
        ## Did a PHEN sub occur?
        phen.sub <- FALSE
        if(phen.node[[edges[e, 1]]] != phen.node[[edges[e, 2]]]) phen.sub <- TRUE
        # SP.dec %in% colnames(Q)[toSub] ## Did EITHER a phen sub AND/OR a snps sub occur?
        
        # print("EDGE"); print(e)
        ## If a sub has occurred on branch e, add it to the n.subs counter
        if(phen.sub){
          N.SUBS.COUNTER <- N.SUBS.COUNTER+1
          # print("TRUE")
        }
        # print("N.SUBS.COUNTER"); print(N.SUBS.COUNTER)
        
      } # end for (e) loop
      
      
      ## STORE FIRST SNPS.ASSOC:
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))
      
      ## STORE PHEN (FOR ALL NODES):
      phen.nodes <- as.vector(unlist(phen.node))
      names(phen.nodes) <- c(1:length(phen.nodes))
      
      phen <- phen.nodes[1:n.ind]
      
      # ## TEMP -- CHECK w PLOT:
      # plot_phen(tree, phen.nodes=phen.nodes, snp.nodes=snps.assoc.nodes[[i]])
      # cor(as.numeric(phen.nodes[1:n.ind]), as.numeric(snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- length(which(phen.nodes[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.OVERLAP[[i]] <- N.overlap
      #######################################################################
      
      
      ## CHECK THAT MIN GRP.SIZE >= THRESHOLD ##
      if(!is.null(grp.min)){
        tab <- table(phen)
        grp.thresh <- (n.ind)*grp.min
        if(min(tab) < grp.thresh){
          toRepeat <- TRUE
        }else{
          toRepeat <- FALSE
        }
      }else{
        toRepeat <- FALSE
      }
      
    } # end WHILE LOOP #########
    
    print("N.SUBS.COUNTER"); print(N.SUBS.COUNTER)
    print("N.overlap"); print(N.overlap[[i]])
    
    ### TEMP -- CHECK:
    # phen.node.ori <- phen.node
    # phen.rec <- as.vector(unlist(phen.node))
    # snp.rec <- as.vector(unlist(snps.assoc.nodes[[3]]))
    # plot_phen(tree, phen.nodes=phen.rec, snp.nodes=snp.rec)
    # title("set3_1 phen vs. snps.assoc3", line=0)
    #
    # cor(as.numeric(snp.rec[1:n.ind]), as.numeric(phen.rec[1:n.ind])) # 0.46 0.63 0.48
    # length(which(as.numeric(snp.rec[1:n.ind]) == as.numeric(phen.rec[1:n.ind])))/n.ind # 0.74 0.81 0.75
    ######
    
    ## get snps.loci for the REMAINING ASSOCIATED snps (ie. 2:n.assoc, conditional on phen) ##
    for(i in 2:n.snps.assoc){
      ## get nt for root at this locus:
      root.nt <- gen.root[snps.assoc[i]]
      
      snp.node <- as.list(rep(NA, length(unique(as.vector(edges)))))
      
      snp.node[[edges[x[1], 1]]] <- root.nt
      
      ## go from last to first edge in edges:
      for(e in x){
        
        probs <- NULL
        
        P <- Qt[[e]]
        
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == FALSE) probs <- P[1,]
        if(snp.node[[edges[e, 1]]] == FALSE & phen.node[[edges[e, 1]]] == TRUE) probs <- P[2,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == FALSE) probs <- P[3,]
        if(snp.node[[edges[e, 1]]] == TRUE & phen.node[[edges[e, 1]]] == TRUE) probs <- P[4,]
        
        ## Now we KNOW, A PRIORI, the phen.node for the DESCENDANT!
        probs.mod <- replace(probs, which(!as.logical(as.numeric(treeWAS::keepLastN(colnames(Q), 1))) == phen.node[[edges[e, 2]]]), 0)
        SP.dec <- sample(colnames(Q), 1, prob = probs.mod)
        
        S.dec <- as.logical(as.numeric(treeWAS::keepFirstN(SP.dec, 1)))
        names(S.dec) <- NULL
        
        snp.node[[edges[e, 2]]] <- S.dec
        
      } # end for (e) loop
      
      ## STORE SNPS.ASSOC (FOR ALL NODES):
      snps.assoc.nodes[[i]] <- as.vector(unlist(snp.node))
      
      ## Get proportion overlap btw phen and snps.assoc.i:
      N.overlap <- length(which(phen.nodes[1:n.ind] == snps.assoc.nodes[[i]][1:n.ind]))
      N.overlap <- max(N.overlap, (n.ind-N.overlap))
      N.OVERLAP[[i]] <- N.overlap
      
    } # end for loop
    
    ## Bind SNPs.ASSOC into matrix:
    snps.assoc.nodes <- do.call("cbind", snps.assoc.nodes)
    
    
    ###################################################
    ## TEMP -- COMPARE PHEN & ALL SNPS.ASSOC w PLOT: ##
    ###################################################
    par(mfrow=c(2,6))
    treeWAS::plot_phen(tree, phen.nodes=phen.nodes, main.title="phen")
    if(n.snps.assoc >= 10){
      ## just plot first 10:
      for(i in 1:5){
        treeWAS::plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                           main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=0, font.main=1)
      }
      treeWAS::plot_phen(tree, phen.nodes=phen.nodes, main.title="phen")
      for(i in 6:10){
        treeWAS::plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                           main.title=paste("snp.assoc", i, sep=" "))
        title(N.OVERLAP[[i]], line=0, font.main=1)
      }
    }else{
      ## plot up to 5:
      if(n.snps.assoc <= 5){
        par(mfrow=c(1,(n.snps.assoc+1)))
        for(i in 1:n.snps.assoc){
          treeWAS::plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                             main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
      }else{
        ## plot btw 5 and 10:
        for(i in 1:5){
          treeWAS::plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                             main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
        treeWAS::plot_phen(tree, phen.nodes=phen.nodes, main.title="phen")
        for(i in 6:n.snps.assoc){
          treeWAS::plot_phen(tree, phen.nodes=snps.assoc.nodes[,i], RTL = TRUE,
                             main.title=paste("snp.assoc", i, sep=" "))
          title(N.OVERLAP[[i]], line=0, font.main=1)
        }
      }
    }
    
    
    par(mfrow=c(1,1)) # end temp panel plot
    
    gc()
    
  } # end of snps.assoc generation
  
  
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###
  
  
  ###########################################
  ## GET COMPLETE SNPS MATRIX ("genomes"): ##
  ###########################################
  
  ## Create snps matrix:
  snps <- temp
  rm(temp)
  
  ## Attach snps.assoc loci to last column:
  if(!is.null(snps.assoc.nodes)){
    snps <- cbind(snps[,c(1:(ncol(snps)-(n.snps.assoc)))], snps.assoc.nodes[1:n.ind,])
  }
  
  
  ## keep only rows containing terminal individuals:
  snps <- snps[1:n.ind, ]
  
  ####   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ###   ####
  
  
  ##############################
  ## PLOTS & TREECONSTRUCTION ##
  ##############################
  tree.reconstructed <- NULL
  if(heatmap == TRUE || reconstruct!=FALSE){
    
    ## CONVERT TO CHARACTER: ##
    dna <- snps
    dna <- replace(dna, which(dna == TRUE), "a")
    dna <- replace(dna, which(dna == "FALSE"), "t")
    
    ## Get DNAbin object:
    dna <- as.DNAbin(dna)
    # rownames(dna) <- c(1:nrow(snps))
    
    
    #############
    ## HEATMAP ##
    #############
    if(heatmap==TRUE){
      heatmap.DNAbin(dna=dna,
                     dist.dna.model=dist.dna.model)
    }
    
    ##########################################
    ## PLOT 2: RECONSTRUCTING THE PHYLOGENY ##
    ##########################################
    if(reconstruct!=FALSE){
      tree.reconstructed <- tree.reconstruct(dna, # [1:n.ind,]
                                             method=reconstruct,
                                             dist.dna.model=dist.dna.model,
                                             plot=TRUE)
    }
    
    ## Remove unnecessary object:
    rm(dna)
    
  } # end heatmap, treeconstruction
  
  
  ##################
  ## CONVERT SNPS ##
  ##################
  
  ## (NO LONGER) CONVERT TO NUMERIC: ##
  ## Convert from logical to binary SNPs (for terminal nodes only):
  # snps <- replace(snps, which(snps == TRUE), 1)
  
  ## Reassort snps.assoc to new columns:
  if(!is.null(snps.assoc)){
    
    ## update snps.assoc to reflect true loci
    gen.size.final <- ncol(snps)
    snps.assoc.loci.ori <- c((gen.size.final-(n.snps.assoc-1)):gen.size.final)
    
    #########################################
    ## RANDOMIZE SNPS.ASSOC LOCI POSITIONS ##
    #########################################
    
    ## Re-enabled snps.assoc loci "randomization" by
    ## just drawing indices and shuffling the columns accordingly...
    ## draw which SNPs will be associated to the phenotype
    if(!is.null(seed)) set.seed(seed)
    if (is.null(ground.truth)) {
      snps.assoc.loci <- sort(sample(c(1:gen.size.final),
                                     n.snps.assoc,
                                     replace=FALSE))
    } else {
      snps.assoc.loci <- sort(ground.truth$nonnulls)
    }
    
    snps.indices <- c(1:gen.size.final)
    snps.ori <- snps
    
    snps.non.assoc <- snps[,c(1:(gen.size.final-n.snps.assoc))]
    snps.assoc <- snps[,snps.assoc.loci.ori]
    snps.new <- matrix(99, nrow=nrow(snps), ncol=gen.size.final)
    snps.new[,snps.indices[-snps.assoc.loci]] <- snps.non.assoc
    snps.new[,snps.assoc.loci] <- snps.assoc
    snps <- snps.new
    snps.assoc <- snps.assoc.loci
    
  } # end snps.assoc randomization
  
  ###############################
  ## Assign row & column names ##
  ###############################
  
  ## assign/generate row.names
  if(!is.null(row.names)){
    ## If row.names have been provided in args list, assign them:
    if(length(row.names) == nrow(snps)){
      ## match tree$tip.label?
      if(!is.null(tree$tip.label)){
        if(all(row.names %in% tree$tip.label) & all(tree$tip.label %in% row.names)){
          ## REORDER to match tree$tip.labs if possible:
          if(!identical(row.names, tree$tip.label)){
            ord <- match(tree$tip.label, row.names)
            row.names <- row.names[ord]
          }
        }
      }
      rownames(snps) <- row.names
    }else{
      if(is.null(rownames(snps))) rownames(snps) <- c(1:nrow(snps))
    }
  }else{
    ## Otherwise, try to assign rownames(snps) to match tree$tip.label:
    if(!is.null(tree$tip.label)){
      if(length(tree$tip.label) == nrow(snps)){
        rownames(snps) <- tree$tip.label
      }else{
        rownames(snps) <- 1:nrow(snps)
        warning("The length of tree$tip.label was not equal to nrow(snps) being simulated;
              rownames(snps) have been set to 1:N and will not match tree$tip.label.")
      }
    }
  }
  
  ## generate column names:
  colnames(snps) <- 1:ncol(snps)
  
  
  #################################################
  ## SIM SET 2 (complementary clade-wise assoc): ##
  #################################################
  sets <- NULL
  if(!is.null(snps.assoc)){
    if(!is.null(set)){
      if(set == 2){
        
        ## Want to divide tree into 2 sets of clades btw 1/3:2/3 and 1/2:1/2
        clades <- tab <- grp.options <- sets.complete <- list()
        
        min.size <- ceiling((n.ind)*(1/3))
        max.size <- floor((n.ind)*(2/3))
        grp1 <- n.ind
        
        ## get 2 sets of clades:
        
        ##########################################
        ## Pick SETS for ANY TREE (pretty sure) ##
        ##########################################
        dec <- grp <- sets.temp <- sets.complete <- list()
        
        inds <- c(1:(n.ind))
        # new.root <- tree$edge[1,1] # initial root
        new.root <- n.ind+1 # initial root
        
        counter <- 0
        #######################################
        ## WHILE LOOP to get size of clades: ##
        #######################################
        while(grp1 < min.size | grp1 > max.size){
          
          ## get all descendants of root node:
          all.dec <- .getDescendants(tree, node=new.root)
          
          ## get all descendants in first 2 major clades:
          dec[[1]] <- .getDescendants(tree, node=all.dec[1])
          dec[[2]] <- .getDescendants(tree, node=all.dec[2])
          
          ## get terminal inds only:
          sets.temp[[1]] <- dec[[1]][which(dec[[1]] %in% inds)]
          sets.temp[[2]] <- dec[[2]][which(dec[[2]] %in% inds)]
          
          grp[[1]] <- length(sets.temp[[1]])
          grp[[2]] <- length(sets.temp[[2]])
          
          max.grp <- which.max(c(grp[[1]], grp[[2]]))
          new.root <- all.dec[max.grp]
          
          set1 <- sets.temp[[max.grp]]
          
          sets <- rep(2, length(inds))
          sets <- replace(sets, set1, 1)
          names(sets) <- tree$tip.label
          
          counter <- counter+1
          
          grp1 <- grp[[max.grp]]
          
        } # end while loop
        
        
        set1 <- names(sets)[which(sets == 1)]
        set2 <- names(sets)[which(sets == 2)]
        ###########
        
        ########################
        ## MODIFY SNPS.ASSOC: ##
        ########################
        snps.assoc.set1 <- 1:round(length(snps.assoc)/2)
        snps.assoc.set2 <- (round(length(snps.assoc)/2)+1):length(snps.assoc)
        
        ## replace set1 snps with 0 at all inds in clade.set1:
        for(e in 1:length(snps.assoc.set1)){
          # snps[which(rownames(snps) %in% set1), snps.assoc[snps.assoc.set1[e]]] <- 0
          snps[which(rownames(snps) %in% set1), snps.assoc[snps.assoc.set1[e]]] <- FALSE
        }
        ## replace set2 snps with 0 at all inds in clade.set2:
        for(e in 1:length(snps.assoc.set2)){
          # snps[which(rownames(snps) %in% set2), snps.assoc[snps.assoc.set2[e]]] <- 0
          snps[which(rownames(snps) %in% set2), snps.assoc[snps.assoc.set2[e]]] <- FALSE
        }
      }
    }
  } # end sim set 2
  
  
  ##################
  ## get RESULTS: ##
  ##################
  out <- list(snps, snps.assoc, tree.reconstructed, sets, phen, phen.nodes, n.mts)
  names(out) <- c("snps", "snps.assoc", "tree.reconstructed", "sets", "phen", "phen.nodes", "n.mts")
  
  return(out)
  
} # end snp.sim.Q


