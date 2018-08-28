# Contains plotting functions to use for identifying clusters of similar gene expression in the data

# load all other functions
source("setup.R")
source("QualityControl.R")
source("CorrectTechnicalNoise.R")
source("Clustering.R")
source("DifferentialGeneExpression.R")

# Create a map of gene expression
Plot_tSNE <- function(q_array, color_ids) {
  cell_ids <- rownames(q_array)
  coordinates <- Rtsne(q_array, is.distance = TRUE, verbose = TRUE, theta = 0.0)
  coordinates <- data.frame(coordinates$Y, row.names = cell_ids)
  coordinates <- cbind.data.frame(coordinates, cluster_ids)
  colnames(coordinates) <- c("X_value", "Y_value", "Color")
  ggplot(data = coordinates, aes(x = X_value, y = Y_value, colour = Color)) +
    geom_point(size = 1) + 
    xlab("tSNE 1") + ylab("tSNE 2") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
}

SC3_consensus_plot <- sc3_plot_consensus(sc3_input, cluster_number)   # saves the SC3 plot so it can be returned
