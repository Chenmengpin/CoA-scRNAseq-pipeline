# produces cell coordinates, does clustering via SC3 and SCENIC

source("setup.R")
source("Preprocessing.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")

# Cluster the cells based on gene expression
SC3clustering <- function(q_array) {
  sce_array <- SingleCellExperiment(assays = list(logcounts = t(q_array)), 
                              rowData = colnames(q_array), colData = rownames(q_array))
  rowData(sce_array)$feature_symbol <- colnames(q_array)
  sc3output <- sc3(sce_array, gene_filter = FALSE, k_estimator = TRUE, kmeans_nstart = 1000)
  assign("sc3_clusters", sc3output, env = .GlobalEnv)
}

# Convert the SCESet output of SC3 into a consensus plot and assign clusters to metadata
SC3Output <- function(sc3_input, sc3_clusters, m_array) {
  SC3_ids <- as.character(sc3_clusters) 
  m_array <- cbind.data.frame(m_array, SC3_ids)
  assign(paste0("SC3_metadata"), m_array, env = .GlobalEnv)
}

# Cluster based on gene regulation
SCENIC_Export <- function(q_array) {
  q_array <- 2^q_array
  q_array <- q_array - 1  # removes the log2 normalization, converting into gene-summarized counts for SCENIC  
  q_array <- t(q_array) # convert to formatting for GENIE3
  write.csv(q_array, file = "SCENIC_clustering/data/GENIE3_import.csv")
}

# imports the output of pySCENIC back into R
SCENIC_Import <- function(result_name) {
  SCENIC_Output <- read.csv("SCENIC_clustering/SCENIC_export.csv", header = TRUE, row.names = 1)
  SCENIC_Output <- t(SCENIC_Output)
  print("SCENIC result loaded")
  SCENIC_result <- AUCell_exploreThresholds(SCENIC_Output, assignCells=TRUE)
  assign(result_name, SCENIC_result, env = .GlobalEnv)
}
# find binary regulon activity based on AUC scores
SCENIC_BinaryArray <- function(SCENIC_result) {
  cellsAssigned <- lapply(SCENIC_result, function(x) x$assignment)
  assignmentTable <- melt(cellsAssigned, value.name="cell")
  print("Binary assignments generated")
  colnames(assignmentTable)[2] <- "geneSet"
  assignmentMat <- table(assignmentTable[, "geneSet"], assignmentTable[, "cell"])
  binary_regulon_array <- as.data.frame.matrix(assignmentMat)
}
