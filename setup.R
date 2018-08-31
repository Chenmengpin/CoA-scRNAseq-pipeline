# Start up everything before running
rm(list = ls())   # Reset the environment, remove all old values before beginning

# Load all packages needed for the pipeline
library(biomaRt)  # for matching Ensembl/GENCODE IDs to gene names
library(tximport)   # does gene-level estimation for C1 data
library(MAST) # does differential gene expression analysis
library(DropletUtils)   # for 10X data preprocessing
library(scran)  # for quality control of single cell data
library(scater)   # for quality control of single cell data
library(DrImpute)   # for imputing probable dropout values
library(dplyr)  # for data manipulation
library(AUCell)   # step 3 of SCENIC, for identifying gene set enrichment
library(Rtsne)  # produces tSNE coordinates for visualization
library(SC3)  # for clustering on gene expression
library(WGCNA)  # finding weighted gene correlation networks
library(dynamicTreeCut)   # for creating improved hierarchical clusters of gene modules
library(ggplot2)  # for making plots
library(tidyr)  # for manipulating dataframes
library(reshape2) # for manipulating dataframes
library(randomcoloR)  # gives consistent color output

# needed for WGCNA to work
options(stringsAsFactors = FALSE)
# enableWGCNAThreads(nThreads = 8)  # massively speeds up WGCNA, does not yet work in RStudio

# needed to allow data frames to remove NaNs
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan)) 
}

# load functions for everything
source("Preprocessing.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")
source("Validation.R")
source("Clustering.R")
source("DifferentialGeneExpression.R")
source("PlotQC.R")
source("PlotValidation.R")
source("PlotClusters.R")
source("PlotGeneExpression.R")
