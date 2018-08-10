# Initialize dependencies and functions from prior steps

source("setup.R")
source("AssembleArrays.R")

# Cell-level quality control
CellQC <- function(q_array, m_array, id, unfiltered_q_array, unfiltered_m_array) {
  library_size <- rowSums(q_array)
  geneset_size <- rowSums(q_array != 0)
  pass_library_qc <- !isOutlier(library_size, nmads = 3, type = "lower", log = TRUE)  # best to use negated versions for metadata import
  pass_gene_qc <- !isOutlier(geneset_size, nmads = 3, type = "lower", log = TRUE)
  m_array <- cbind.data.frame(m_array, library_size, geneset_size)   # add all library and geneset info to metadata
  q_array <- q_array[pass_library_qc == TRUE & pass_gene_qc == TRUE,]   # dual filters on library and geneset size, needs to pass both
  m_array <- m_array[pass_library_qc == TRUE & pass_gene_qc == TRUE,]
  assign(paste0("cellqc_quant_",id), q_array, env = .GlobalEnv)   # these need to be made like this so that it returns both with custom names
  assign(paste0("GCqc_metadata_",id), m_array, env = .GlobalEnv)
  print("Beginning repeat for metadata annotation")
  library_size <- rowSums(unfiltered_q_array)  # repeat to produce full metadata list containing all cells
  geneset_size <- rowSums(unfiltered_q_array != 0)
  pass_library_qc <- !isOutlier(library_size, nmads = 3, type = "lower", log = TRUE)  
  pass_gene_qc <- !isOutlier(geneset_size, nmads = 3, type = "lower", log = TRUE)
  unfiltered_m_array <- cbind.data.frame(unfiltered_m_array, library_size, pass_library_qc, geneset_size, pass_gene_qc)
  assign(paste0("unfiltered_metadata_",id), unfiltered_m_array, env = .GlobalEnv)
}

# Gene-level quality control
GeneQC <- function(q_array, id) {
  cells_per_gene <- as.vector(colSums(q_array != 0))
  q_array <- q_array[, cells_per_gene >= 3]
  assign(paste0("geneqc_quant_",id), q_array, env = .GlobalEnv)
}

# Scaling by size factor
NormalizeCountData <- function(q_array, m_array, id, unfiltered_q_array, unfiltered_m_array) {
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
  print("Beginning repeat for metadata annotation")   # repeat to produce full metadata list containing all cells
  unfiltered_q_array <- jitter(t(unfiltered_q_array))   # need to jitter to prevent zero rank variances
  deconvolution_clusters <- quickCluster(unfiltered_q_array)
  if (method_id == "10X") {
    pass_size_qc <- computeSumFactors(q_array, sizes = seq(20, 120, 2), # computes size factors per cell, use more pools for higher precision
                                      clusters = deconvolution_clusters, min.mean = .1)   # UMI counts need lower threshold 
  } else {
    pass_size_qc <- computeSumFactors(q_array, sizes = seq(20, 120, 2), # computes size factors per cell, use more pools for higher precision
                                      clusters = deconvolution_clusters)   # clusters improve performance by reducing differential expression 
  }
  print("Size factors computed for metadata")
  pass_size_qc <- pass_size_qc > 0  # creates vector of whether or not the size factor is positive
  unfiltered_m_array <- cbind.data.frame(unfiltered_m_array, pass_size_qc)
  assign(paste0("unfiltered_metadata_",id), unfiltered_m_array, env = .GlobalEnv)
}

# Compensating for batch effects and inter-individual variability
UnifyDatasets <- function(q_array, m_array, unfiltered_m_array, dataset_names) {
  cell_ids <- lapply(q_array, rownames)   # these ensure cell and gene info survives matrix transformation
  names(cell_ids) <- paste0(dataset_names,"_names")   # each element needs a name before export to outside environment
  names(q_array) <- dataset_names   # each element needs a name before export to outside environment
  gene_ids <- lapply(q_array, colnames)
  commongenes <- Reduce(intersect, gene_ids)  # need to create list of all genes present
  for (i in 1:length(q_array)) {
    q_array[[i]] <- q_array[[i]][1:nrow(q_array[[i]]), intersect(colnames(q_array[[i]]), commongenes)]
    print(i)
  }
  q_array <- lapply(q_array, t)   # transforms dataframes into a group of matrices compatible with mixed nearest neighbor correction
  list2env(q_array, environment())  # splits the list into a group of matrices
  print("Finished prep for MNN clustering")
  q_array <- mnnCorrect(get(dataset_names))   # correct for batch effects WAIT AND FIX THIS 
  print("Finished MNN clustering")
  q_array <- as.data.frame(q_array$corrected)   # unify the batch corrected quant arrays
  m_array <- bind_rows(m_array)   # unify metadata into single arrays
  cell_ids <- rownames(m_array)
  unfiltered_m_array <- bind_rows(unfiltered_m_array)   # rip the row IDs from the metadata, because rows are changed from list conversion in quant array
  assign(unified_metadata, m_array, env = .GlobalEnv)
  assign(unified_unfiltered_metadata, unfiltered_m_array, env = .GlobalEnv)
  q_array <- as.data.frame(t(q_array), row.names = cell_ids)  # make the quant array into a dataframe and label it for compatibility
  colnames(q_array) <- commongenes
  assign(unified_quants, q_array, env = .GlobalEnv)
}

ImputeDropouts <- function(q_array, id) {
  cell_ids <- rownames(q_array)   # these ensure cell and gene info survives matrix transformation
  gene_ids <- colnames(q_array)   
  q_array <- t(q_array)   # DrImpute requires a matrix with cells in columns and genes in rows
  q_array <- DrImpute(q_array, ks = 10:30, # maximize number of k-values in case dataset is highly heterogeneous
                      dists = c("spearman", "pearson", "euclidean")) # maximize number of distributions for better consensus
  q_array <- data.frame(t(q_array), row.names = cell_ids) # returns to metadata-compatible format
  colnames(q_array) <- gene_ids
  assign(paste0("imputed_quant_",id), q_array, env = .GlobalEnv)  # returns original quant array identifier with modifier indicating normalization
}