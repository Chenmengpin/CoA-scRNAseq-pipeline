# Start up everything before running
rm(list = ls())   # Reset the environment, remove all old values before beginning

# Load all packages needed for the pipeline
library(biomaRt)  # for matching Ensembl/GENCODE IDs to gene names
library(dplyr)  # for the pipe operator
