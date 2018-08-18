# produces cell coordinates, does clustering via SC3 and SCENIC

source("setup.R")
source("Preprocessing.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")

# Create initial cell coordinates based on gene expression
tSNECoordinates <- function(q_array, cell_ids) {
  coordinates <- Rtsne(q_array, is.distance = TRUE, verbose = TRUE, theta = 0.0)
  coordinates <- data.frame(coordinates$Y, row.names = cell_ids)
  colnames(coordinates) <- c("X", "Y")
  assign("tSNECoordinates", coordinates, env = .GlobalEnv)
}

# Cluster the cells based on gene expression
SC3clustering <- function(q_array) {
  sce_array <- SingleCellExperiment(assays = list(logcounts = t(q_array)), 
                              rowData = colnames(q_array), colData = rownames(q_array))
  rowData(sce_array)$feature_symbol <- colnames(q_array)
  sc3output <- sc3(sce_array, gene_filter = FALSE, k_estimator = TRUE, kmeans_nstart = 1000)
  assign("sc3_clusters", sc3output, env = .GlobalEnv)
}

# Convert the SCESet output of SC3 into a consensus plot and assign clusters to metadata
SC3Output <- function(sc3_clusters, clusters, m_array) {
  SC3_ids <- as.character(sc3_input) 
  m_array <- cbind.data.frame(m_array, SC3_ids)
  assign(paste0("SC3_metadata"), m_array, env = .GlobalEnv)
  SC3_consensus_plot <- sc3_plot_consensus(sc3_input, 22)   # saves the SC3 plot so it can be returned
  return(SC3_consensus_plot)
}

# Cluster based on gene regulation
SCENICGeneInference <- function(q_array) {
  q_array <- 2^q_array
  q_array <- q_array - 1  # removes the log2 normalization, converting into gene-summarized counts for SCENIC  
  q_array <- t(q_array)
  gn_array <- GENIE3(q_array, verbose = TRUE)
}
