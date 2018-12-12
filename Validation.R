# produces list of validation metrics in one step

# run the validation script
ValidateQC <- function(q_array, m_array) {
  q_array <- q_array[is.na(m_array$Batch) == FALSE,]  # removes any possible NA values that could ruin analysis
  m_array <- m_array[is.na(m_array$Batch) == FALSE,]
  dim_number <- sum(length(unique(m_array$Batch)), choose(length(unique(m_array$Batch)), 2), 
                    length(unique(m_array$Sample)), choose(length(unique(m_array$Sample)), 2),
                    length(unique(m_array$Region)), choose(length(unique(m_array$Region)), 2),
                    length(unique(m_array$Method)), choose(length(unique(m_array$Method)), 2))  # tallies up total number of types of pairwise comparisons to make
  ValidationCorrelation <- data.frame(matrix(0, nrow = 100, ncol = dim_number * 6))   # creates target matrix to put values into
  ValidationGenes <- data.frame(matrix(0, nrow = 100, ncol = dim_number * 6))
  iteration_number <- 1
  for (type in 1:4) {
    type_id <- colnames(m_array)[type]
    for (i in 1:length(unique(m_array[, type]))) {
      i_label <- unique(m_array[, type])[i]
      for (j in 1:length(unique(m_array[, type]))) {
        j_label <- unique(m_array[, type])[j]
        if (is.na(match(paste(type_id, j_label, "v", i_label, "1v1"), colnames(ValidationCorrelation))) == TRUE) {
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 1)  # identify samples to be randomly selected
            cell_1_match <- which(rownames(m_array) == cell_1)
            cell_1_array <- unlist(q_array[cell_1_match,])
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 1)
            cell_2_match <- which(rownames(m_array) == cell_2)
            cell_2_array <- unlist(q_array[cell_2_match,])
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- as.numeric(cell_1_array + cell_2_array)
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "1v1", i_label, j_label, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v1")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v1")
          iteration_number <- iteration_number + 1
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 1)  # identify samples to be randomly selected
            cell_1_match <- which(rownames(m_array) == cell_1)
            cell_1_array <- unlist(q_array[cell_1_match,])
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 10)
            cell_2_match <- match(cell_2, rownames(m_array))
            cell_2_array <- q_array[cell_2_match,]
            cell_2_sum <- as.numeric(colSums(cell_2_array))
            cell_2_array <- as.numeric(colMeans(cell_2_array))
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- cell_1_array + cell_2_sum
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "1v10", i_label, j_label, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v10")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v10")
          iteration_number <- iteration_number + 1
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 1)  # identify samples to be randomly selected
            cell_1_match <- which(rownames(m_array) == cell_1)
            cell_1_array <- unlist(q_array[cell_1_match,])
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 100)
            cell_2_match <- match(cell_2, rownames(m_array))
            cell_2_array <- q_array[cell_2_match,]
            cell_2_sum <- as.numeric(colSums(cell_2_array))
            cell_2_array <- as.numeric(colMeans(cell_2_array))
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- cell_1_array + cell_2_sum
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "1v100", i_label, j_label, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v100")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "1v100")
          iteration_number <- iteration_number + 1
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 10)
            cell_1_match <- match(cell_1, rownames(m_array))
            cell_1_array <- q_array[cell_1_match,]
            cell_1_sum <- as.numeric(colSums(cell_1_array))
            cell_1_array <- as.numeric(colMeans(cell_1_array))
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 10)
            cell_2_match <- match(cell_2, rownames(m_array))
            cell_2_array <- q_array[cell_2_match,]
            cell_2_sum <- as.numeric(colSums(cell_2_array))
            cell_2_array <- as.numeric(colMeans(cell_2_array))
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- cell_1_sum + cell_2_sum
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "10v10", i_label, j_label, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "10v10")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "10v10")
          iteration_number <- iteration_number + 1
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 10)
            cell_1_match <- match(cell_1, rownames(m_array))
            cell_1_array <- q_array[cell_1_match,]
            cell_1_sum <- as.numeric(colSums(cell_1_array))
            cell_1_array <- as.numeric(colMeans(cell_1_array))
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 100)
            cell_2_match <- match(cell_2, rownames(m_array))
            cell_2_array <- q_array[cell_2_match,]
            cell_2_sum <- as.numeric(colSums(cell_2_array))
            cell_2_array <- as.numeric(colMeans(cell_2_array))
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- cell_1_sum + cell_2_sum
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "10v100", i_label, j_label, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "10v100")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "10v100")
          iteration_number <- iteration_number + 1
          for (k in 1:100) {
            cell_1 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[i]], size = 100)
            cell_1_match <- match(cell_1, rownames(m_array))
            cell_1_array <- q_array[cell_1_match,]
            cell_1_sum <- as.numeric(colSums(cell_1_array))
            cell_1_array <- as.numeric(colMeans(cell_1_array))
            cell_2 <- sample(x = rownames(m_array)[m_array[, type] == unique(m_array[, type])[j]], size = 100)
            cell_2_match <- match(cell_2, rownames(m_array))
            cell_2_array <- q_array[cell_2_match,]
            cell_2_sum <- as.numeric(colSums(cell_2_array))
            cell_2_array <- as.numeric(colMeans(cell_2_array))
            ValidationCorrelation[k, iteration_number] <- cor(x = cell_1_array, y = cell_2_array)
            sum_vector <- cell_1_sum + cell_2_sum
            ValidationGenes[k, iteration_number] <- length(sum_vector[sum_vector > .02])
            print(paste(type_id, "100v100", i, j, k))
          }
          colnames(ValidationCorrelation)[iteration_number] <- paste(type_id, i_label, "v", j_label, "100v100")
          colnames(ValidationGenes)[iteration_number] <- paste(type_id, i_label, "v", j_label, "100v100")
          iteration_number <- iteration_number + 1
        } else { next }
      }
    } 
  }
  ValidationCorPVal <- data.frame(matrix(0, ncol = length(colnames(ValidationCorrelation)), nrow = length(colnames(ValidationCorrelation)),
                                         dimnames = list(colnames(ValidationCorrelation), colnames(ValidationCorrelation))))
  ValidationGenePVal <- data.frame(matrix(0, ncol = length(colnames(ValidationGenes)), nrow = length(colnames(ValidationGenes)),
                                          dimnames = list(colnames(ValidationCorrelation), colnames(ValidationCorrelation))))  
  for (l in 1:length(colnames(ValidationCorPVal))) {
    for (m in 1:length(colnames(ValidationCorPVal))) {
      Cor_Test <- t.test(x = ValidationCorrelation[, l], ValidationCorrelation[, m], paired = FALSE)
      ValidationCorPVal[l,m] <- Cor_Test$p.value
      Gene_Test <- t.test(x = ValidationCorrelation[, l], ValidationCorrelation[, m], paired = FALSE)
      ValidationGenePVal[l,m] <- Gene_Test$p.value
      print(paste(l, m))
    }
    ValidationCorPVal[,m] <- p.adjust(ValidationCorPVal[,m], method = "bonferroni")
    ValidationGenePVal[,m] <- p.adjust(ValidationGenePVal[,m], method = "bonferroni")
  }
  ValidationQC <- list(ValidationCorrelation, ValidationGenes, ValidationCorPVal, ValidationGenePVal)
  names(ValidationQC) <- c("Correlations", "Genes", "CorPVal", "GenePVal")
  assign("ValidationQC", ValidationQC, env = .GlobalEnv)
}
