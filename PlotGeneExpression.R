# Create plots of all summarized gene expression and its interaction with all clusters

# load all other functions
source("setup.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")
source("Clustering.R")
source("DifferentialGeneExpression.R")

# Quality control for soft thresholding
# 2 is SFT R-squared, 5 is mean connectivity, 8 is density, 9 is centralization, 10 is heterogeneity
WGCNA_Plot_SoftThreshold <- function(thresholds, metric_id, metric_label) {
  thresholds_data <- thresholds$fitIndices
  thresholds_data$Color <- thresholds_data$Power == thresholds$powerEstimate
  image <- ggplot(data = thresholds_data, aes(x = thresholds_data$Power, y = thresholds_data[, metric_id], colour = thresholds_data$Color)) +
    geom_line() + geom_point(size = 2) +
    scale_colour_manual(values = c("black", "red")) +
    scale_x_continuous(expand = c(.02, 0)) + scale_y_continuous(expand = c(.02, 0)) +
    xlab("Soft Threshold (Power)") + ylab(metric_label) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="none")
  return(image)
}

WGCNA_Plot_ModuleTree <- function(module_tree, modules) {
  image <- plotDendroAndColors(module_tree, modules$colors, "Modules", dendroLabels = FALSE,
                      hang = 0.03, addGuide = TRUE, guideHang = 0.05)
  dev.off()
  return(image)
}

# plots the relationship between WGCNA gene modules and clusters in terms of R-squared
WGCNA_Plot_ModuleTraitSignificance <- function(modules, cluster_ids) {
  module_names <- modules$dendro$labels
  cluster_names <- levels(cluster_ids)
  cluster_ids <- as.numeric(cluster_ids)
  cluster_matrix <- matrix(0, nrow = length(cluster_ids), ncol = max(cluster_ids))
  cluster_matrix[col(cluster_matrix) == cluster_ids] <- 1
  ModulePopSignificance <- corAndPvalue(x = modules$newMEs, y = cluster_matrix, alternative = "two.sided")
  image <- labeledHeatmap(ModulePopSignificance$cor, colorLabels = FALSE, colors = matlab.like(50),
                          xLabels = cluster_names, yLabels = rownames(ModulePopSignificance$cor), 
                          xLabelsAngle = 90, setStdMargins = FALSE)  
  dev.off()
  return(image)
}

# plots the relationship between WGCNA gene modules with one another in terms of R-squared
WGCNA_Plot_ModuleTraitSignificance <- function(modules) {
  module_names <- modules$dendro$labels
  ModulePopSignificance <- corAndPvalue(x = modules$newMEs, y = modules$newMEs, alternative = "two.sided")
  image <- labeledHeatmap(ModulePopSignificance$cor, colorLabels = FALSE, colors = matlab.like(50),
                          xLabels = rownames(ModulePopSignificance$cor), 
                          yLabels = rownames(ModulePopSignificance$cor), 
                          xLabelsAngle = 90, setStdMargins = FALSE)  
  dev.off()
  return(image)
}

# plots relationship between WGCNA gene modules in terms of hierarchical clustering
WGCNA_Plot_ModuleHierarchy <- function(modules) {
  module_tree <- modules$dendro
  image <- plot(modules_tree, xlab = "", sub = "")
  dev.off()
  return(image)
}


