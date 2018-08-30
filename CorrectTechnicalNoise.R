# Remove technical noise in each dataset

# Compensating for batch effects and inter-individual variability
PrepDatasets <- function(..., dataset_names, m_array_list, qc_m_array_list, gene_m_array_list) {
  q_array <- list(...)
  cell_ids <- lapply(q_array, rownames)   # these ensure cell and gene info survives matrix transformation
  combined_cell_ids <<- Reduce(c, cell_ids) # concatenates all cell_id vectors in the list together
  names(cell_ids) <- paste0(dataset_names,"_names")   # each element needs a name before export to outside environment
  names(q_array) <- dataset_names   # each element needs a name before export to outside environment
  gene_ids <- lapply(q_array, colnames)
  commongenes <<- Reduce(intersect, gene_ids)  # need to create list of all genes present, save for batch correction step
  for (i in 1:length(q_array)) { # subsets each quant array by the set of genes all have in common
    q_array[[i]] <- q_array[[i]][1:nrow(q_array[[i]]), intersect(colnames(q_array[[i]]), commongenes)]
    print(i)
  }
  q_array <- lapply(q_array, t)   # transforms dataframes into a group of matrices compatible with mixed nearest neighbor correction
  for (i in 1:length(q_array)) { # gives each dataframe its name and and names the genes still present
    rownames(q_array[[i]]) <- commongenes 
    assign(dataset_names[i], q_array[[i]], envir = .GlobalEnv)
    print(i)
  }
  m_array <- bind_rows(m_array_list)   # unify metadata into single arrays
  rownames(m_array) <- combined_cell_ids
  qc_m_array <- bind_rows(qc_m_array_list)   # rip the row IDs from the metadata, because rows are changed from list conversion in quant array
  cells_per_gene <- Reduce("+", lapply(gene_m_array_list, "[[", 3))   # finds the sum of cells across all batches expressing the gene 
  pass_cellnumber_qc_every_batch <- Reduce(intersect, lapply(gene_m_array_list, function(x) x$pass_cellnumber_qc))
  for (i in 1:length(gene_m_array_list)) {  # gets total expression of each gene in each sample
    gene_m_array_list[[i]][1] <- gene_m_array_list[[i]][1] * gene_m_array_list[[i]][3]
    all_genes <- rownames(gene_m_array_list[[i]])
    print(i)
  }
  mean_nontransformed_expression <- Reduce("+", lapply(gene_m_array_list, "[[", 1))   # finds the sum of expression across all cells in all batches
  mean_nontransformed_expression <- mean_nontransformed_expression / cells_per_gene   # gets mean expression across all batches
  mean_transformed_expression <- log2(mean_nontransformed_expression + 1)   # transforms all mean expression patterns
  gene_ids <- lapply(gene_m_array_list, function(x) x[,4]) # turns all pass_cellnumber_qc vectors into their own list
  pass_cellnumber_qc_every_batch <- apply(data.frame(gene_ids), 1, all) # matches the true and false values in each vector, all need to be TRUE to mark TRUE
  gene_m_array <- cbind.data.frame(mean_nontransformed_expression, mean_transformed_expression, cells_per_gene, pass_cellnumber_qc_every_batch)
  rownames(gene_m_array) <- all_genes
  assign("unified_metadata", m_array, env = .GlobalEnv)
  assign("unified_qc_metadata", qc_m_array, env = .GlobalEnv)
  assign("unified_gene_metadata", gene_m_array, env = .GlobalEnv)
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
  assign(paste0("imputed_quants",id), q_array, env = .GlobalEnv)  # returns original quant array identifier with modifier indicating normalization
}

# Filter down to highly variable genes for SC3 and edgeR
FilterHVG <- function(q_array, spike_ins, gene_m_array) {
  q_array_cv2 <- t(q_array)
  hvg_array <- improvedCV2(q_array_cv2, is.spike = spike_ins, log.prior = 1)
  q_array <- q_array[, hvg_array$p.value < .05]
  assign("HVG_quants", q_array, env = .GlobalEnv)
  gene_m_array$row <- rownames(gene_m_array)    # add a 'helper' column for the merge function.
  hvg_array$row <- rownames(hvg_array)
  gene_m_array <- merge(gene_m_array, hvg_array, by='row', all=TRUE)
  gene_m_array <- gene_m_array[, -c(1, 6, 7, 9)]
  assign("HVGs_unified_gene_metadata", gene_m_array, env = .GlobalEnv)
}
