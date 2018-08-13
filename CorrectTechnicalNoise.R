# Initialize dependencies and functions from prior steps

source("setup.R")
source("Preprocessing.R")
source("QualityControl.R")

# Compensating for batch effects and inter-individual variability
PrepDatasets <- function(..., dataset_names, m_array, qc_m_array) {
  q_array <- list(...)
  cell_ids <- lapply(q_array, rownames)   # these ensure cell and gene info survives matrix transformation
  combined_cell_ids <<- Reduce(c, cell_ids)
  names(cell_ids) <- paste0(dataset_names,"_names")   # each element needs a name before export to outside environment
  names(q_array) <- dataset_names   # each element needs a name before export to outside environment
  gene_ids <- lapply(q_array, colnames)
  commongenes <<- Reduce(intersect, gene_ids)  # need to create list of all genes present, save for batch correction step
  for (i in 1:length(q_array)) {
    q_array[[i]] <- q_array[[i]][1:nrow(q_array[[i]]), intersect(colnames(q_array[[i]]), commongenes)]
    print(i)
  }
  q_array <- lapply(q_array, t)   # transforms dataframes into a group of matrices compatible with mixed nearest neighbor correction
  for (i in 1:length(q_array)) { # need to 
    rownames(q_array[[i]]) <- commongenes 
    assign(dataset_names[i], q_array[[i]], envir = .GlobalEnv)
    print(i)
  }
  m_array <- bind_rows(m_array)   # unify metadata into single arrays
  rownames(m_array) <- combined_cell_ids
  qc_m_array <- bind_rows(qc_m_array)   # rip the row IDs from the metadata, because rows are changed from list conversion in quant array
  assign("unified_metadata", m_array, env = .GlobalEnv)
  assign("unified_qc_metadata", qc_m_array, env = .GlobalEnv)
}

BatchCorrect <- function(..., commongenes, combined_cell_ids) {
  q_array <- mnnCorrect(...)
  q_array <- as.data.frame(q_array$corrected)   # unify the batch corrected quant arrays
  q_array <- as.data.frame(t(q_array), row.names = combined_cell_ids)  # make the quant array into a dataframe and label it for compatibility
  colnames(q_array) <- commongenes
  assign("unified_quants", q_array, env = .GlobalEnv)
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