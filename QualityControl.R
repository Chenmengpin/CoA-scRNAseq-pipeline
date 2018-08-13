# Initialize dependencies and functions from prior steps

source("setup.R")
source("Preprocessing.R")

# Cell-level quality control
CellQC <- function(q_array, m_array, id, qc_m_array) {
  library_size <- rowSums(q_array)
  geneset_size <- rowSums(q_array != 0)
  library_qc <- !isOutlier(library_size, nmads = 3, type = "lower", log = TRUE)  # best to use negated versions for metadata import
  gene_qc <- !isOutlier(geneset_size, nmads = 3, type = "lower", log = TRUE)
  m_array <- cbind.data.frame(m_array, library_size, geneset_size)   # add all library and geneset info to metadata
  q_array <- q_array[library_qc == TRUE & gene_qc == TRUE,]   # dual filters on library and geneset size, needs to pass both
  library_m_array <- m_array[library_qc == TRUE,]
  m_array <- library_m_array[gene_qc == TRUE,]
  assign(paste0("cellqc_quant_",id), q_array, env = .GlobalEnv)   # these need to be made like this so that it returns both with custom names
  assign(paste0("GCqc_metadata_",id), m_array, env = .GlobalEnv)
  print("Beginning metadata QC annotation")
  pass_library_qc <- is.element(rownames(qc_m_array), rownames(library_m_array))  
  pass_gene_qc <- is.element(rownames(qc_m_array), rownames(m_array))
  qc_m_array <- cbind.data.frame(qc_m_array, pass_library_qc, pass_gene_qc)
  assign(paste0("QC_metadata_",id), qc_m_array, env = .GlobalEnv)
}

# Gene-level quality control
GeneQC <- function(q_array, id) {
  cells_per_gene <- as.vector(colSums(q_array != 0))
  q_array <- q_array[, cells_per_gene >= 3]
  assign(paste0("geneqc_quant_",id), q_array, env = .GlobalEnv)
}

# Scaling by size factor
NormalizeCountData <- function(q_array, m_array, id, qc_m_array) {
  method_id <- m_array[1,4]
  cell_ids <- rownames(q_array)   # these ensure cell and gene info survives matrix transformation
  gene_ids <- colnames(q_array)   
  q_array <- t(q_array)   # scran normalization requires a matrix with cells in columns and genes in rows
  deconvolution_clusters <- quickCluster(q_array)   # low-level clustering improves deconvolution performance by minimizing differential expression
  print("Deconvolution clusters calculated")
  if (method_id == "10X") {
    size_factors <- computeSumFactors(q_array, sizes = seq(20, 120, 2), # computes size factors per cell, use more pools for higher precision
                                      clusters = deconvolution_clusters, min.mean = .1)   # UMI counts need lower threshold 
  } else {
    size_factors <- computeSumFactors(q_array, sizes = seq(20, 120, 2), # computes size factors per cell, use more pools for higher precision
                                      clusters = deconvolution_clusters)   # clusters improve performance by reducing differential expression 
  }
  print("Size factors computed for QC")
  q_array <- scale(q_array, center = FALSE, scale = size_factors)   # performs scaling with these factors
  print("Cells scaled by library size")
  q_array <- data.frame(t(q_array), row.names = cell_ids) # returns to metadata-compatible format
  colnames(q_array) <- gene_ids
  q_array <- q_array[size_factors > 0,]
  m_array <- m_array[size_factors > 0,]
  q_array <- q_array + 1  # prevents undefined values for zeroes in log transformation
  q_array <- log2(q_array)  # log-transforms data to account for heteroscedasticity, log2 used because it is fine-grained and easy to represent fold changes
  assign(paste0("normalized_quant_",id), q_array, env = .GlobalEnv)  # returns original quant array identifier with modifier indicating normalization
  assign(paste0("normalized_metadata_",id), m_array, env = .GlobalEnv)
  print("Beginning metadata QC annotation")
  pass_size_qc <- is.element(rownames(qc_m_array), rownames(m_array))
  qc_m_array <- cbind.data.frame(qc_m_array, pass_size_qc)
  assign(paste0("QC_metadata_",id), qc_m_array, env = .GlobalEnv)
}