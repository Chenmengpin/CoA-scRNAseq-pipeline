# Start up everything before running
rm(list = ls())   # Reset the environment, remove all old values before beginning

# Load all packages needed for the pipeline
library(biomaRt)  # for matching Ensembl/GENCODE IDs to gene names
library(tximport)   # does gene-level estimation for C1 data
library(DropletUtils)   # for 10X data preprocessing
library(scran)  # for quality control of single cell data
library(scater)   # for quality control of single cell data
library(DrImpute)   # for imputing probable dropout values
library(dplyr)  # for data manipulation
library(GENIE3) # step 1 of SCENIC, does gene network inference
library(RcisTarget)   # step 2 of SCENIC, identifies transcription factor binding motifs in genes
library(AUCell)   # step 3 of SCENIC, for identifying gene set enrichment
library(Rtsne)  # produces tSNE coordinates for visualization
library(SC3)  # for clustering on gene expression
