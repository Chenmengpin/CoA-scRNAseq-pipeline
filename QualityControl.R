# Initialize dependencies and functions from prior steps

source("setup.R")
source("AssembleArrays.R")

# WORK IN PROGRESS, MORE NEEDS TO BE DONE HERE #

# TXIMPORT STEP FOR C1 DATA

# Cell-level quality control
CellQC <- function(q_array, m_array, q_id, m_id) {
  library_size <- rowSums(q_array)
  geneset_size <- rowSums(q_array != 0)
  pass_library_qc <- !isOutlier(library_size, nmads = 3, type = "lower", log = TRUE)  # best to use negated versions for metadata import
  pass_gene_qc <- !isOutlier(geneset_size, nmads = 3, type = "lower", log = TRUE)
  m_array <- cbind.data.frame(m_array, library_size, pass_library_qc, geneset_size, pass_gene_qc)   # add all library and geneset info to metadata
  m_array_all <- m_array
  q_array <- q_array[pass_library_qc == TRUE & pass_gene_qc == TRUE,]   # dual filters on library and geneset size, needs to pass both
  m_array <- m_array[pass_library_qc == TRUE & pass_gene_qc == TRUE,]
  assign(paste0("cellqc_",q_id), q_array, env = .GlobalEnv)   # these need to be made like this so that it returns both with custom names
  assign(paste0("qc_",m_id), m_array, env = .GlobalEnv)
  assign(paste0("unfiltered_",m_id), m_array_all, env = .GlobalEnv)
}

# Gene-level quality control
GeneQC <- function(q_array, q_id) {
  cells_per_gene <- as.vector(colSums(q_array != 0))
  q_array <- q_array[, cells_per_gene >= 3]
  assign(paste0("geneqc_",q_id), q_array, env = .GlobalEnv)
}

# Scaling by size factor
NormalizeCountData <- function(q_array, m_array, q_id) {
  batch_id <- m_array[1,1]
  sample_id <- m_array[1,2]
  method_id <- m_array[1,3]
  cell_ids <- rownames(q_array)   # these ensure cell and gene info survives matrix transformation
  gene_ids <- colnames(q_array)   
  q_array <- t(q_array)   # scran normalization requires a matrix with cells in columns and genes in rows
  print("Quant array transformed to matrix")
  deconvolution_clusters <- quickCluster(q_array)   # low-level clustering improves deconvolution performance by minimizing differential expression
  print("Deconvolution clusters calculated")
  if (method_id == "10X") {
    size_factors <- computeSumFactors(q_array, sizes = seq(20, 120, 2), # computes size factors per cell, use more pools for higher precision
                                      min.mean = .1)   # UMI counts need lower threshold 
  } else {
    size_factors <- computeSumFactors(q_array, sizes = seq(20, 120, 2) # computes size factors per cell, use more pools for higher precision
                                      )   # clusters improve performance by reducing differential expression 
  }
  print("Size factors computed")
  q_array <- scale(q_array, center = FALSE, scale = size_factors)   # performs scaling with these factors
  print("Cells scaled by library size")
  q_array <- data.frame(t(q_array), row.names = cell_ids) # returns to metadata-compatible format
  colnames(q_array) <- gene_ids
  for (i in 1:nrow(q_array)) {  # use rows to reduce number of loops
    q_array[i,] <- q_array[i,] + 1  # prevents undefined values for zeroes in log transformation
    q_array[i,] <- log2(q_array[i,])   # log-transforms data to account for heteroscedasticity, log2 used because it is fine-grained and easy to represent fold changes
    print(i)
  }
  assign(paste0("normalized_",q_id), q_array, env = .GlobalEnv)  # returns original quant array identifier with modifier indicating normalization
}