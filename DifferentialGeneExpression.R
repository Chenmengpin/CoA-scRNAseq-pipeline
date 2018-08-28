# Figures out differential gene expression between all different clusters

# load all other functions
source("setup.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")
source("Clustering.R")

# does differential gene expression analysis
DGE_edgeR <- function(q_array, cluster_id) {
  DEG_PValue <- list()
  DEG_logFC <- list()
  q_array <- 2^q_array
  q_array <- q_array - 1  # removes the log2 normalization, converting into gene-summarized counts for edgeR
  q_array <- t(q_array)
  q_array[q_array < 0] <- 0   # if there are negative values, it will crash. all negatives are near-zero values from normalization
  DGE_array <- DGEList(counts = q_array, group = cluster_id)
  DGE_array <- calcNormFactors(DGE_array)
  print("Beginning model fitting")
  cell_det_rate <- scale(colMeans(q_array > 0))
  genewise_models <- model.matrix(~ cell_det_rate + cluster_id)
  model <- estimateDisp(DGE_array, design = genewise_models, robust = TRUE)
  fit_model <- glmQLFit(model, design = genewise_models, robust = TRUE)
  print("Model fitting done")
  for (i in 1:nlevels(cluster_id)) {
    DEresults_logFC <- data.frame(rownames(q_array))
    rownames(DEresults_logFC) <- DEresults_logFC[,1]
    DEresults_PValue <- data.frame(rownames(q_array))
    rownames(DEresults_PValue) <- DEresults_PValue[,1]
    DEresults_coln <- "Gene"
    print(paste("Set of comparisons for cluster", i))
    selected_cluster <- as.character(i)
    for (j in 1:nlevels(cluster_id)) {
      print(paste("Comparison against cluster", j))
      if (j == selected_cluster) { next }
      contrast <- numeric(ncol(genewise_models))
      contrast[i] <- 1
      contrast[j] <- -1
      print("Running tests")
      qlf_test <- glmQLFTest(fit_model, contrast = contrast)
      print("Tests done")
      newCompPV <- as.vector(qlf_test$table$PValue)
      DEresults_PValue <- cbind.data.frame(DEresults_PValue, newCompPV)
      newCompFC <- as.vector(qlf_test$table$logFC)
      DEresults_logFC <- cbind.data.frame(DEresults_logFC, newCompFC)
      DEresults_coln <- c(DEresults_coln, paste0(i, "_vs_", j))
    }
    colnames(DEresults_PValue) <- DEresults_coln
    DEresults_PValue <- DEresults_PValue[,-1] 
    colnames(DEresults_logFC) <- DEresults_coln
    DEresults_logFC <- DEresults_logFC[,-1] 
    for (m in 1:nrow(DEresults_PValue)) {
      print(m)
      DEresults_PValue[m,] <- p.adjust(DEresults_PValue[m,], method = "BH")
    }
    DEresults_PValue$SigComps <- rowSums(DEresults_PValue < .05)  
    DEG_PValue[[i]] <- DEresults_PValue
    names(DEG_PValue)[i] <- paste0("PValue_Cluster_", i)
    DEG_logFC[[i]] <- DEresults_logFC
    names(DEG_logFC)[i] <- paste0("logFC_Cluster_", i)
  }
  assign("DEG_PValue", DEG_PValue, env = .GlobalEnv)
  assign("DEG_logFC", DEG_logFC, env = .GlobalEnv)
  assign("DEG_model", model, env = .GlobalEnv)
  assign("DEG_QLfitting", fit_model, env = .GlobalEnv)
}

# find variable gene co-expression modules
DGE_WGCNA <- function(q_array) {
  threshold_list <- pickSoftThreshold(q_array, powerVector = seq(1, 30, by = 1), 
                                      corFnc = bicor, networkType = "signed hybrid", moreNetworkConcepts = TRUE, verbose = 5)
  assign("WGCNA_thresholds", file = threshold_list, env = .GlobalEnv)
  adjacency_matrix <- adjacency(q_array, type = "signed hybrid", power = threshold_list$powerEstimate, corFnc = bicor)
  TOMatrix <- TOMsimilarity(adjacency_matrix, TOMType = "signed", verbose = 5)
  inverse_TOMatrix <- 1 - TOMatrix
  module_tree <- hclust(as.dist(inverse_TOMatrix), method = "average")
  module_ids <- cutreeDynamic(dendro = module_tree, distM = inverse_TOMatrix, deepSplit = 4)
  module_ids <- mergeCloseModules(exprData = q_array, colors = as.numeric(module_ids), corFnc = bicor, verbose = 5)
  assign("WGCNA_hierarchy", file = module_tree, env = .GlobalEnv)
  assign("WGCNA_gene_clusters", file = module_ids, env = .GlobalEnv)
}

# test modules identified in WGCNA for reproducibility IN PROGRESS
test_WGCNA <- function() {
  
}
