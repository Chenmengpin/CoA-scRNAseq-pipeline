# Contains plotting functions to make all plots for use in validation figures

# plot gene expression correlation across samples
# 1 for correlation, 2 for number of genes
Plot_Validation_Boxplots <- function(ValidationList, metric, label_list) {
  metric_id <- match(metric, names(ValidationList))
  plot_dataset <- ValidationList[[metric_id]]
  if (metric_id == 1) {
    y_label <- "R-value"
  } else {
    y_label <- "Genes"
  }
  plot_dataset <- melt(plot_dataset)
  group_type <- rep(1:(length(rownames(plot_dataset))/6), each = 600, length.out = length(rownames(plot_dataset)))
  group_type <- as.character(group_type)
  analysis_type <- rep(1:6, each = 100, length.out = length(rownames(plot_dataset)))
  analysis_type <- as.character(analysis_type)
  plot_dataset <- cbind.data.frame(plot_dataset, analysis_type, group_type)
  tiff("yes.tiff", width = 9000, height = 9000, res = 300)
  ggplot(plot_dataset, aes(x = group_type, y= value, colour = analysis_type)) + geom_point(position = position_jitterdodge(jitter.width = 0.1)) + geom_boxplot() +
    scale_colour_manual(name = "Comparison", labels = c("1 v 1", "1 v 10", "1 v 100", "10 v 10", "10 v 100", "100 v 100"), values = c("blue", "red", "forestgreen", "purple", "orange", "pink")) +
    xlab("") + ylab(y_label) + scale_x_discrete(labels = label_list) + guides(colour = guide_legend(nrow = 1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.position = "bottom", legend.key=element_blank(), legend.box = "horizontal",
          axis.line = element_line(colour = "black"))
  dev.off()
}

# plot the representative graphs for each validation metric
Plot_Validation_Scatter <- function(q_array, m_array, type, specific_x, specific_y, cell_num_x, cell_num_y) {
  type_id <- match(type, colnames(m_array))
  q_array <- q_array[is.na(m_array$Batch) == FALSE,]  # removes any possible NA values that could ruin analysis
  m_array <- m_array[is.na(m_array$Batch) == FALSE,]
  if (cell_num_x == 1) {
    cell_x <- sample(x = rownames(m_array)[m_array[, type_id] == specific_x], size = 1)  # identify samples to be randomly selected
    cell_x_match <- which(rownames(m_array) == cell_x)
    cell_x_array <- unlist(q_array[cell_x_match,])
    cell_x_label <- "Cell, Log(Counts + 1)"
  } else {
    cell_x <- sample(x = rownames(m_array)[m_array[, type_id] == specific_x], size = cell_num_x)
    cell_x_match <- match(cell_x, rownames(m_array))
    cell_x_array <- q_array[cell_x_match,]
    cell_x_array <- as.numeric(colMeans(cell_x_array))
    cell_x_label <- "Cells, Mean(Log(Counts + 1))"
  }
  if (cell_num_y == 1) {
    cell_y <- sample(x = rownames(m_array)[m_array[, type_id] == specific_y], size = 1)  # identify samples to be randomly selected
    cell_y_match <- which(rownames(m_array) == cell_y)
    cell_y_array <- unlist(q_array[cell_y_match,]) 
    cell_y_label <- "Cell, Log(Counts + 1)"
  } else {
    cell_y <- sample(x = rownames(m_array)[m_array[, type_id] == specific_y], size = cell_num_y)
    cell_y_match <- match(cell_y, rownames(m_array))
    cell_y_array <- q_array[cell_y_match,]
    cell_y_array <- as.numeric(colMeans(cell_y_array))
    cell_y_label <- "Cells, Mean(Log(Counts + 1))"
  }
  plot_array <- cbind.data.frame(cell_x_array, cell_y_array)
  r_value <- round(cor(cell_x_array, cell_y_array), digits = 2)
  ggplot(plot_array, aes(x = cell_x_array, y = cell_y_array)) + 
    geom_point(size = 1) + geom_smooth(method = "glm", color = "red", fullrange = TRUE) + 
    geom_text(aes(x = min(cell_x_array), y = max(cell_y_array)), label = paste("R =", r_value), vjust = "inward", hjust = "inward") +
    xlab(paste(type, specific_x, "-", cell_num_x, cell_y_label)) + 
    ylab(paste(type, specific_y, "-", cell_num_y, cell_y_label)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
}

# plot the p-values in a heatmap for each validation metric
Plot_Validation_Heatmap <- function(ValidationList, metric) {
  metric_id <- match(metric, names(ValidationList))
  plot_dataset <- ValidationList[[metric_id]]
  if (metric_id == 3) {
    y_label <- "R-value"
  } else {
    y_label <- "Genes"
  }
  plot_dataset <- melt(plot_dataset)
  group_type <- as.character(rep(unique(plot_dataset$variable), times = nlevels(plot_dataset$variable)))
  plot_dataset <- cbind.data.frame(plot_dataset, group_type) 
  group_type <- as.character(rep(unique(plot_dataset$variable), each = nlevels(plot_dataset$variable)))
  plot_dataset$variable <- group_type
  ggplot(plot_dataset, aes(x= variable, y = group_type, fill = -log10(value))) + 
    geom_tile() + xlab("") + ylab("") + 
    scale_fill_gradientn(colours = matlab.like(10), name = "-log10(adj. p-value)") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.ticks = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top")
}