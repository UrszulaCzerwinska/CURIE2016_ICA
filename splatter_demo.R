## ----knitr-options, echo = FALSE, message = FALSE, warning = FALSE-------
# To render an HTML version that works nicely with github and web pages, do:
# rmarkdown::render("vignettes/splatter.Rmd", "all")
knitr::opts_chunk$set(fig.align = 'center', fig.width = 6, fig.height = 5,
                      dev = 'png')

## ----install, eval = FALSE-----------------------------------------------
#  source("https://bioconductor.org/biocLite.R")
#  biocLite("splatter")

## ----install-github, eval = FALSE----------------------------------------
#  biocLite("Oshlack/splatter", dependencies = TRUE,
#           build_vignettes = TRUE)

## ----quickstart----------------------------------------------------------
# Load package
library(splatter)

# Load example data
data("sc_example_counts")
# Estimate parameters from example data
params <- splatEstimate(sc_example_counts)
# Simulate data using estimated parameters
sim <- splatSimulate(params, dropout.present = FALSE)

## ----SplatParams---------------------------------------------------------
params <- newSplatParams()
params

## ----getParam------------------------------------------------------------
getParam(params, "nGenes")

## ----setParam------------------------------------------------------------
params <- setParam(params, "nGenes", 5000)
getParam(params, "nGenes")

## ----getParams-setParams-------------------------------------------------
# Set multiple parameters at once (using a list)
params <- setParams(params, update = list(nGenes = 8000, mean.rate = 0.5))
# Extract multiple parameters as a list
getParams(params, c("nGenes", "mean.rate", "mean.shape"))
# Set multiple parameters at once (using additional arguments)
params <- setParams(params, mean.shape = 0.5, de.prob = 0.2)
params

## ----newSplatParams-set--------------------------------------------------
params <- newSplatParams(lib.loc = 12, lib.scale = 0.6)
getParams(params, c("lib.loc", "lib.scale"))

## ----splatEstimate-------------------------------------------------------
# Check that sc_example counts is an integer matrix
class(sc_example_counts)
typeof(sc_example_counts)
# Check the dimensions, each row is a gene, each column is a cell
dim(sc_example_counts)
# Show the first few entries
sc_example_counts[1:5, 1:5]

params <- splatEstimate(sc_example_counts)

## ----splatSimulate-------------------------------------------------------
sim <- splatSimulate(params, nGenes = 1000, dropout.present = FALSE)
sim

## ----SCESet--------------------------------------------------------------
# Access the counts
counts(sim)[1:5, 1:5]
# Information about genes
head(fData(sim))
# Information about cells
head(pData(sim))
# Gene by cell matrices
names(assayData(sim))
# Example of cell means matrix
get_exprs(sim, "CellMeans")[1:5, 1:5]

## ----pca-----------------------------------------------------------------
plotPCA(sim)

## ----groups--------------------------------------------------------------
sim.groups <- splatSimulate(groupCells = c(100, 100), method = "groups",
                            verbose = FALSE)
plotPCA(sim.groups, colour_by = "Group")

## ----paths---------------------------------------------------------------
sim.paths <- splatSimulate(method = "paths", verbose = FALSE)
plotPCA(sim.paths, colour_by = "Step")

## ----listSims------------------------------------------------------------
listSims()

## ----listSims-table------------------------------------------------------
knitr::kable(listSims(print = FALSE))

## ----lengths-------------------------------------------------------------
sim <- simpleSimulate(verbose = FALSE)
sim <- addGeneLengths(sim)
head(fData(sim))

## ----TPM-----------------------------------------------------------------
tpm(sim) <- calculateTPM(sim, fData(sim)$Length)
tpm(sim)[1:5, 1:5]

## ----comparison----------------------------------------------------------
sim1 <- splatSimulate(nGenes = 1000, groupCells = 20, verbose = FALSE)
sim2 <- simpleSimulate(nGenes = 1000, nCells = 20, verbose = FALSE)
comparison <- compareSCESets(list(Splat = sim1, Simple = sim2))

names(comparison)
names(comparison$Plots)

## ----comparison-means----------------------------------------------------
comparison$Plots$Means

## ----comparison-libsize-features-----------------------------------------
library("ggplot2")
ggplot(comparison$PhenoData,
       aes(x = total_counts, y = total_features, colour = Dataset)) +
    geom_point()

## ----sessionInfo---------------------------------------------------------
sessionInfo()

