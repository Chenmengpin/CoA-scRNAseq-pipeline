# Figures out differential gene expression between all different clusters


# does differential gene expression analysis via MAST
DGE_MAST <- function(q_array, cluster_id, gene_metadata) {
  MAST_clusters <- as.character(c(cluster_id, 0))   # sequester intercept away from real data, where we have no logFC for this
  intercept_vector <- rep(0, times = length(colnames(q_array))) # treat intercept vector like a cell with its own cluster
  MAST_array <- rbind.data.frame(q_array, intercept_vector)
  MAST_DGE_list <- list()
  MAST_GO_export <- list()
  MAST_array <- t(MAST_array)   # need this to format for analysis via MAST
  cell_det_rate <- scale(colMeans(MAST_array > 0))  # CDR depends on finding the mean number of genes detected per cell
  MAST_sca <- FromMatrix(exprsArray = MAST_array, 
                         fData = data.frame(primerid = rownames(MAST_array)), # primerid is what MAST calls the names of genes
                         cData = data.frame(wellKey = colnames(MAST_array), # wellKey is what MAST calls cells
                                            cluster = MAST_clusters, cdr = cell_det_rate))  # these are the co-variates of interest
  zlm_result <- zlm(~cdr + cluster, sca = MAST_sca)   # need to compute the model for MAST
  genes_all <- rownames(gene_metadata)
  GO_frame <- rep(1, length = length(genes_all))  # creating an array for import into topGO
  names(GO_frame) <- genes_all
  for (i in 1:(length(colnames(zlm_result@coefC)) - 2)) {   # need to do this to prevent namelength errors
    contrast_value <- colnames(zlm_result@coefC)[i+2]   # first two slots filled by the CDR and intercept, which are not going to be considered
    LRT_array <- summary(zlm_result, doLRT = contrast_value)    # this calculates the p-values for each
    LRT_dt <- LRT_array$datatable   # removes the superficial list layer
    LRT_dt <- LRT_dt[LRT_dt$contrast == contrast_value]   # limits analysis to the cluster of interest
    P_values <- LRT_dt[LRT_dt$component == "H",]    # use the p-value for the hurdle model
    P_values$`Pr(>Chisq)` <- p.adjust(P_values$`Pr(>Chisq)`, method = "BH")   # need to correct P-values, else way too likely to find signficace
    gene_list <- P_values$primerid
    pval_list <- P_values$`Pr(>Chisq)`
    LFC <- LRT_dt[LRT_dt$component == "logFC",]   # limits to the log-fold change
    LFC_list <- LFC$coef
    filtered_array <- cbind.data.frame(gene_list, pval_list, LFC_list)  # this is the array we actually care about
    colnames(filtered_array) <- c("Gene", "Pvalue", "logFC")
    filtered_array$logFC[is.nan(filtered_array$logFC)] <- 0
    filtered_array$PvalSig <- filtered_array$Pvalue < .01   # cuts off significant p-value at FDR-corrected 0.01
    filtered_array$logFCSig <- filtered_array$logFC > 0   # sets baseline for minimum increase to be considered a marker gene
    filtered_array$marker <- filtered_array$PvalSig & filtered_array$logFCSig == TRUE   # unifies marker gene criteria
    MAST_DGE_list[[i]] <- filtered_array
    cluster_GO_frame <- GO_frame
    cluster_GO_genes <- filtered_array$Gene[filtered_array$marker == TRUE]
    cluster_GO_frame[names(cluster_GO_frame) %in% cluster_GO_genes] <- 2
    cluster_GO_frame <- factor(cluster_GO_frame, levels = c(1, 2), labels = c("0", "1"))
    MAST_GO_export[[i]] <- cluster_GO_frame
    print(i)
  }
  names(MAST_DGE_list) <- colnames(zlm_result@coefC)[3:length(colnames(zlm_result@coefC))]
  names(MAST_GO_export) <- names(MAST_DGE_list)
  assign("MAST_DGE_array", MAST_DGE_list, env = .GlobalEnv)
  assign("MAST_GO_export", MAST_GO_export, env = .GlobalEnv)
}

# find variable gene co-expression modules
DGE_WGCNA <- function(q_array, gene_metadata) {
  threshold_list <- pickSoftThreshold(q_array, powerVector = seq(1, 30, by = 1), 
                                      corFnc = bicor, networkType = "signed hybrid", moreNetworkConcepts = TRUE, verbose = 5)
  adjacency_matrix <- adjacency(q_array, type = "signed hybrid", power = threshold_list$powerEstimate, corFnc = bicor)
  TOMatrix <- TOMsimilarity(adjacency_matrix, TOMType = "signed", verbose = 5)
  inverse_TOMatrix <- 1 - TOMatrix
  print("Refining module delineation")
  module_tree <- hclust(as.dist(inverse_TOMatrix), method = "average")
  module_ids <- cutreeDynamic(dendro = module_tree, distM = inverse_TOMatrix, deepSplit = 4)
  module_ids_all <- mergeCloseModules(exprData = q_array, colors = as.numeric(module_ids), corFnc = bicor, verbose = 5)
  module_ids <- module_ids_all$colors
  print("Creating WGCNA module gene array for GO export")
  gene_array <- colnames(q_array)
  WGCNA_GO_export <- list()
  genes_all <- rownames(gene_metadata)
  GO_frame <- rep(0, length = length(genes_all))  # creating an array for import into topGO
  names(GO_frame) <- genes_all
  for (i in 1:max(module_ids)) {
    module_GO_frame <- GO_frame
    module_GO_genes <- gene_array[module_ids == i]
    module_GO_frame[names(module_GO_frame) %in% module_GO_genes] <- 1
    module_GO_frame <- factor(module_GO_frame)
    WGCNA_GO_export[[i]] <- module_GO_frame
    print(i)
  }
  names(WGCNA_GO_export) <- paste0("module", 1:max(module_ids))
  #CytoScape_WGCNA <- exportNetworkToCytoscape(adjacency_matrix, nodeAttr = module_array$Module,
  #                                            edgeFile = "WGCNA_edge", nodeFile = "WGCNA_node")
  assign("WGCNA_thresholds", file = threshold_list, env = .GlobalEnv)
  assign("WGCNA_hierarchy", file = module_tree, env = .GlobalEnv)
  assign("WGCNA_gene_clusters", file = module_ids, env = .GlobalEnv)
  assign("WGCNA_GO", file = 'WGCNA_GO_export', env = .GlobalEnv)
}
