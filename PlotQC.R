# Contains plotting functions to make all plots for use in QC figures

# plot amount of cells eliminated by each quality control measure
PlotCellExclusionQC <- function(qc_m_array, group_id, colorlist) {
  qc_m_array <- na.omit(qc_m_array)
  summary <- qc_m_array %>%
                group_by(qc_m_array[, group_id]) %>%
                summarize(Initial <- length(Sample), 
                          EmptyDroplets <- sum(pass_emptydrop_qc),
                          GenesetSize <- sum(pass_gene_qc),
                          LibrarySize <- sum(pass_library_qc),
                          MitoFraction <- sum(mt_qc),
                          SizeFactors <- sum(pass_size_qc),
                          Final <- sum(pass_size_qc))
  colnames(summary) <- c("Group", "0", "1", "2", "3", "4", "5", "6")
  summary <- melt(summary)
  summary$variable <- as.numeric(summary$variable) - 1
  ggplot(data = summary, aes(x = variable, y = value, colour = Group)) + 
    geom_step(na.rm = TRUE, size = 1) + 
    scale_colour_manual(values = colorlist, name = colnames(qc_m_array[group_id])) +
    scale_x_continuous(expand = c(0, 0), labels = c("Initial", "Empty Droplets", "Gene Detection", "Library Size", "Mitochondrial Fraction", "Size Factors", "Final")) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, (max(summary$value) * 1.1))) +
    xlab("Quality control checks") + ylab("Number of cells passing QC") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"))
}

# plot library sizes that were excluded based on QC
PlotUMICountQC <- function(library_array) {
  library_array$library_qc <- ifelse(library_array$library_qc == TRUE, "Pass", "Low MAD outlier")
  if (length(unique(library_export$library_qc)) == 1) {
    ggplot(data = library_array, aes(library_size, fill = library_qc)) +
      geom_histogram(binwidth = 500, colour = "black") +
      scale_fill_manual(values = c("darkgrey", "red"), name = "UMI Count QC Check") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      xlab("Total UMI count (non-normalized)") + ylab("Number of cells") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom") 
  }
  else {
    ggplot(data = library_array, aes(library_size, fill = library_qc)) +
      geom_histogram(binwidth = 500, colour = "black") +
      scale_fill_manual(values = c("red", "darkgrey"), name = "UMI Count QC Check") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      xlab("Total UMI count (non-normalized)") + ylab("Number of cells") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom") 
  }
}

# plot geneset sizes that were excluded based on QC
PlotGenesetQC <- function(geneset_array) {
  geneset_array$gene_qc <- ifelse(geneset_array$gene_qc == TRUE, "Pass", "Less than 1500 genes")
  if (length(unique(library_export$library_qc)) == 1) {
    ggplot(data = geneset_array, aes(geneset_size, fill = gene_qc)) +
      geom_histogram(binwidth = 50, colour = "black") +
      scale_fill_manual(values = c("darkgrey", "red"), name = "Gene Detection QC Check") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      xlab("Number of genes detected") + ylab("Number of cells") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom")
  }
  else {
    ggplot(data = geneset_array, aes(geneset_size, fill = gene_qc)) +
      geom_histogram(binwidth = 50, colour = "black") +
      scale_fill_manual(values = c("red", "darkgrey"), name = "Gene Detection QC Check") +
      scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
      xlab("Number of genes detected") + ylab("Number of cells") + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
            legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom") 
  }
}

PlotMitoQC <- function(mt_array) {
  mt_array$mt_qc <- ifelse(mt_array$mt_qc == TRUE, "Pass", "More than 12% mtRNA")
  ggplot(data = mt_array, aes(mt_fraction, fill = mt_qc)) +
    geom_histogram(binwidth = .01, colour = "black") +
    scale_fill_manual(values = c("red", "darkgrey"), name = "Mitochondrial QC Check") +
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) +
    xlab("Proportion of counts derived from mtRNA") + ylab("Number of cells") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom")
}

# plot comparison between size of library and size factor
PlotLibrarySizeFactorQC <- function(qc_m_array) {
  qc_m_array$size_factor <- qc_m_array$size_factor + .01
  qc_m_array$pass_size_qc <- ifelse(qc_m_array$pass_size_qc == TRUE, "Pass", "Non-positive size factor")
  ggplot(data = qc_m_array, aes(x = size_factor, y = library_size, colour = pass_size_qc)) + 
    geom_point(na.rm = TRUE, size = 1) + 
    geom_point(aes(x = mean(size_factor), y = mean(library_size)), colour = "black", size = 3) +
    scale_colour_manual(values = c("red", "darkgrey"), name = "Size Factor QC Check") +
    scale_x_log10(expand = c(.01, 0)) + scale_y_log10(expand = c(.01, 0)) +
    xlab("Size factor (X + 1e-2)") + ylab("Total read depth") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"), legend.position="bottom")
}

# plot comparison between size of library and geneset size
PlotLibraryGenesQC <- function(m_array, group_id, colorlist) {
  m_array <- na.omit(m_array)
  summary <- m_array %>%
                group_by(m_array[, group_id]) %>%
                summarize(Mean_Library = mean(library_size),
                          Mean_Geneset = mean(geneset_size))
  colnames(summary)[1] <- "Group"
  ggplot(data = m_array, aes(x = library_size, y = geneset_size, colour = m_array[, group_id])) + 
    geom_point(na.rm = TRUE, size = 1) + 
    geom_point(data = summary, aes(x = Mean_Library, y = Mean_Geneset, fill = Group), colour = "black", pch = 21, size = 5) +
    scale_colour_manual(values = colorlist, name = colnames(m_array[group_id])) +
    scale_fill_manual(values = colorlist, name = colnames(m_array[group_id])) +
    scale_x_continuous(expand = c(.01, 0)) + scale_y_continuous(expand = c(.01, 0)) +
    xlab("Total UMI count (non-normalized)") + ylab("Genes detected") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"))
}

# plot amount of genes passing each hurdle for each successive analysis
PlotGeneExclusionQC <- function(BV_gene_array) {
  starting_genes <- length(BV_gene_array$mean_nontransformed_expression)
  cells_per_gene <- sum(BV_gene_array$cells_per_gene >= 3)
  common_genes <- sum(BV_gene_array$pass_cellnumber_qc_every_batch)
  BV_gene_array[is.na(BV_gene_array)] <- 1
  bio_variable <- sum(BV_gene_array$pass_bv_qc)
  final_genes <- bio_variable
  points <- c(0, 1, 2, 3, 4)
  nums <- c(starting_genes, cells_per_gene, common_genes, bio_variable, final_genes)
  summary <- cbind.data.frame(points, nums)
  ggplot(data = summary, aes(x = points, y = nums)) + 
    geom_step(na.rm = TRUE, size = 1) + 
    scale_x_continuous(expand = c(0, 0), limits = c(0, (max(points) * 1.03)), labels = c("Initial", "At Least 3 Cells", "In All Datasets", "Biologically Variable", "Final")) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, (max(summary$nums) * 1.1))) +
    xlab("Quality control checks") + ylab("Number of genes passing QC") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"))
}

# plot the relationship between cells expressing a gene and its average expression
PlotExpressionCellRelationship <- function(HVG_gene_array) {
  HVG_gene_array[is.nan(HVG_gene_array)] <- 0
  colnames(HVG_gene_array)[4] <- "Sufficient"
  HVG_gene_array$Sufficient <- ifelse(HVG_gene_array$Sufficient == TRUE, "Yes", "No")
  ggplot(data = HVG_gene_array, aes(x = mean_transformed_expression, y = cells_per_gene, colour = Sufficient)) +
    geom_point(size = 1) +
    #geom_smooth(aes(group = identity), method = "loess", colour = "red", se = TRUE) +
    scale_colour_manual(values = c("red", "black")) +
    scale_x_log10(expand = c(.01, 0)) + scale_y_log10(expand = c(.01, 0)) +
    xlab("Mean normalized read count") + ylab("Number of cells") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"))
}

# plot the genes that are highly variable
PlotHVGs <- function(HVG_gene_array) {
  HVG_gene_array <- na.omit(HVG_gene_array)
  colnames(HVG_gene_array)[3:4] <- c("Significance", "identity")
  HVG_gene_array$Significance <- HVG_gene_array$p.value < .05
  HVG_gene_array$Significance <- ifelse(HVG_gene_array$Significance == TRUE, "p < 0.05", "p > 0.05")
  ggplot(data = HVG_gene_array, aes(x = mean_transformed_expression, y = cv2, colour = Significance)) +
    geom_point(size = 1) +
    geom_smooth(aes(group = identity), method = "loess", colour = "red", se = TRUE) +
    scale_colour_manual(values = c("red", "black")) +
    scale_x_log10(expand = c(.01, 0)) + scale_y_log10(expand = c(.01, 0)) +
    xlab("Mean normalized read count") + ylab("Squared coefficient of variation") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          legend.key=element_blank(), axis.line = element_line(colour = "black"))
}
