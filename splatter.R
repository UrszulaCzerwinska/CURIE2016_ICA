source("https://bioconductor.org/biocLite.R")
biocLite("scater")
#install.packages("checkmate")
library(splatter)

###following tutorial
data("sc_example_counts")
dim(sc_example_counts)
summary(sc_example_counts)
head(sc_example_counts)
# Estimate parameters from example data
params <- splatEstimate(sc_example_counts)
# Simulate data using estimated parameters
sim <- splatSimulate(params, dropout.present = FALSE)
head(sim)
sim.data <- counts(sim)

head(sim.data )
dim(sim.data )
plotPCA(sim)

#compare if covariance matrix is preserved
cov.original <- cov(t(sc_example_counts))
cov.simul <- cov(t(sim.data ))

library("RSpectra")
#eigen of both
eigen.original <-eigs_sym(cov.original, 100)
eigen.simul <-eigs_sym(cov.simul,100)

#plot 
barplot(eigen.original$values, main= "eigenvalues oroginal")
barplot(eigen.simul$values,  main= "eigenvalues simul")
plot(eigen.original$vectors[,1], eigen.simul$vectors[,1], pch=16) #idk how to interpret this 
plot(eigen.original$vectors[,2], eigen.simul$vectors[,2], pch=16) #idk how to interpret this 


###check with my data

summary(unlog_data(compute_trans_cell_type(caf)))
caf.counts <- as.matrix(unlog_data(compute_trans_cell_type(caf)))
params <- splatEstimate(caf.counts)
sim.caf <- splatSimulate(params, dropout.present = FALSE)
head(sim)

# Get the parameters we are going to use
nCells <- getParam(params, "nCells")
nGenes <- getParam(params, "nGenes")
nGroups <- getParam(params, "nGroups")
group.cells <- getParam(params, "groupCells")


method <- "single"
cell.names <- paste0("Cell", seq_len(nCells))
gene.names <- paste0("Gene", seq_len(nGenes))


dummy.counts <- matrix(1, ncol = nCells, nrow = nGenes)
rownames(dummy.counts) <- gene.names
colnames(dummy.counts) <- cell.names
phenos <- new("AnnotatedDataFrame", data = data.frame(Cell = cell.names))
rownames(phenos) <- cell.names
features <- new("AnnotatedDataFrame", data = data.frame(Gene = gene.names))
rownames(features) <- gene.names
sim <- newSCESet(countData = dummy.counts, phenoData = phenos,
                 featureData = features)

splatSimLibSizes <- function(sim, params) {
  
  nCells <- getParam(params, "nCells")
  lib.loc <- getParam(params, "lib.loc")
  lib.scale <- getParam(params, "lib.scale")
  
  exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
  pData(sim)$ExpLibSize <- exp.lib.sizes
  
  return(sim)
}
sim <- splatSimLibSizes(sim, params)

getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {
  
  is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
  n.selected <- sum(is.selected)
  dir.selected <- (-1) ^ rbinom(n.selected, 1, neg.prob)
  facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
  # Reverse directions for factors that are less than one
  dir.selected[facs.selected < 1 & dir.selected == -1] <- 1
  dir.selected[facs.selected < 1 & dir.selected == 1] <- -1
  factors <- rep(1, n.facs)
  factors[is.selected] <- facs.selected ^ dir.selected
  
  return(factors)
}
splatSimGeneMeans <- function(sim, params) {

sim=sim
params=params
    
  nGenes <- getParam(params, "nGenes")
  mean.shape <- getParam(params, "mean.shape")
  mean.rate <- getParam(params, "mean.rate")
  out.prob <- getParam(params, "out.prob")
  out.facLoc <- getParam(params, "out.facLoc")
  out.facScale <- getParam(params, "out.facScale")
  
  # Simulate base gene means
  base.means.gene <- rgamma(nGenes, shape = mean.shape, rate = mean.rate)
  
  # Add expression outliers
  outlier.facs <- getLNormFactors(nGenes, out.prob, 0, out.facLoc,
                                  out.facScale)
  median.means.gene <- median(base.means.gene)
  outlier.means <- median.means.gene * outlier.facs
  is.outlier <- outlier.facs != 1
  means.gene <- base.means.gene
  means.gene[is.outlier] <- outlier.means[is.outlier]
  
  fData(sim)$BaseGeneMean <- base.means.gene
  fData(sim)$OutlierFactor <- outlier.facs
  fData(sim)$GeneMean <- means.gene
  
  return(sim)
}

sim <- splatSimGeneMeans(sim, params)
splatSimSingleCellMeans <- function(sim, params) {
  
  nCells <- getParam(params, "nCells")
  cell.names <- pData(sim)$Cell
  gene.names <- fData(sim)$Gene
  exp.lib.sizes <- pData(sim)$ExpLibSize
  
  cell.means.gene <- as.matrix(fData(sim)[, rep("GeneMean", nCells)])
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  set_exprs(sim, "BaseCellMeans") <- base.means.cell
  
  return(sim)
}
sim <- splatSimSingleCellMeans(sim, params)

summary(round(getParam(params, "mean.rate"),2))

splatSimBCVMeans <- function(sim, params) {
  
  cell.names <- pData(sim)$Cell
  gene.names <- fData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  bcv.common <- getParam(params, "bcv.common")
  bcv.df <- getParam(params, "bcv.df")
  base.means.cell <- get_exprs(sim, "BaseCellMeans")
  
  bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
    sqrt(bcv.df / rchisq(nGenes, df = bcv.df)) #creates Inf
  
  means.cell <- matrix(rgamma(nGenes * nCells, shape = 1 / (bcv ^ 2),
                              scale = base.means.cell * (bcv ^ 2)),
                       nrow = nGenes, ncol = nCells) #Warning message: In rgamma(nGenes * nCells, shape = 1/(bcv^2), scale = base.means.cell *  : NAs produced
  
  colnames(means.cell) <- cell.names
  rownames(means.cell) <- gene.names
  
  set_exprs(sim, "BCV") <- bcv
  set_exprs(sim, "CellMeans") <- means.cell
  
  return(sim)
}

sim <- splatSimBCVMeans(sim, params)

splatSimTrueCounts <- function(sim, params) {
  
  cell.names <- pData(sim)$Cell
  gene.names <- fData(sim)$Gene
  nGenes <- getParam(params, "nGenes")
  nCells <- getParam(params, "nCells")
  cell.means <- get_exprs(sim, "CellMeans")
  
  true.counts <- matrix(rpois(nGenes * nCells, lambda = cell.means),
                        nrow = nGenes, ncol = nCells)
  
  colnames(true.counts) <- cell.names
  rownames(true.counts) <- gene.names
  
  set_exprs(sim, "TrueCounts") <- true.counts
  
  return(sim)
}
sim <- splatSimTrueCounts(sim, params)
sim <- splatSimDropout(sim, params)

sce <- newSCESet(countData = counts(sim),
                 phenoData = new("AnnotatedDataFrame", data = pData(sim)),
                 featureData = new("AnnotatedDataFrame", data = fData(sim)))

# Add intermediate matrices stored in assayData
for (assay.name in names(assayData(sim))) {
  if (!(assay.name %in% names(assayData(sce)))) {
    set_exprs(sce, assay.name) <- get_exprs(sim, assay.name)
  }
}


splatSimSingleCellMeans <- function(sim, params) {
  
  nCells <- getParam(params, "nCells")
  cell.names <- pData(sim)$Cell
  gene.names <- fData(sim)$Gene
  exp.lib.sizes <- pData(sim)$ExpLibSize
  
  cell.means.gene <- as.matrix(fData(sim)[, rep("GeneMean", nCells)])
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  set_exprs(sim, "BaseCellMeans") <- base.means.cell
  
  return(sim)
}



#########
summary(unlog_data(compute_trans_cell_type(t_cells)))
t_cells.counts <- as.matrix(unlog_data(compute_trans_cell_type(t_cells)))
params_tc <- splatEstimate(t_cells.counts)
sim.t_cells <- splatSimulate(params_tc, dropout.present = FALSE)
head(sim.t_cells)