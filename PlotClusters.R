# Contains plotting functions to use for identifying clusters of similar gene expression in the data

# Create a map of gene expression
Plot_tSNE <- function(q_array, color_ids) {
  umap_config <- umap.defaults
  umap_config$n.neighbors <- 30
  umap_config$n.epochs <- 1000
  umap_config$input <- 'dist'
  umap_config$verbose <- TRUE
  coordinates <- umap(q_array)
  coordinates <- data.frame(coordinates$Y, row.names = cell_ids)
  coordinates <- cbind.data.frame(coordinates, color_ids)
  colnames(coordinates) <- c("X_value", "Y_value", "Color")
  ggplot(data = coordinates, aes(x = X_value, y = Y_value, colour = Color)) +
    geom_point(size = 1) + 
    xlab("tSNE 1") + ylab("tSNE 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

SC3_consensus_plot <- sc3_plot_consensus(sc3_input, cluster_number)   # saves the SC3 plots so it can be returned
SC3_marker_genes <- sc3_plot_markers(sc3_clusters, k = 9) ####FIX THIS

Plot_BinaryArray <- function(binary_regulon_array, cluster_id) {
  binary_regulon_heatmap <- binary_regulon_array
  binary_regulon_heatmap <- data.frame(t(binary_regulon_array))
  binary_regulon_heatmap$cluster <- cluster_id
  binary_regulon_heatmap <- binary_regulon_heatmap[order(binary_regulon_heatmap$cluster),]
  cluster_ann <- data.frame(binary_regulon_heatmap$cluster)
  binary_regulon_heatmap <- binary_regulon_heatmap[, - grep("cluster", colnames(binary_regulon_heatmap))]
  binary_regulon_heatmap <- t(binary_regulon_heatmap)
  color_list <- list(cluster_ann = distinctColorPalette(k = nlevels(cluster_id)))
  tiff("yep.tiff", width = 900, height = 900, res = 300)
  aheatmap(binary_regulon_heatmap, annCol = cluster_ann, annColors = color_list,
           Colv = NA, scale = "none", color = "black", legend = FALSE, annLegend = TRUE)
  dev.off()
}
