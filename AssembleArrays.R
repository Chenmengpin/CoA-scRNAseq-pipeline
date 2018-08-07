# Matches cell IDs and annotations to values and assigns gene identities

source("setup.R")   # Import dependencies

# Create the quant array for gene expression by each cell
AssembleQuantArray <- function(quant_csv, gene_ids, batch_id, sample_id) {
  transcript_ids <- read.table(gene_ids, colClasses = "character")  # import the annotations from alignment
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")   # connect to Ensembl via BioMart
  ann2gene <- getBM(c("ensembl_gene_id_version", "external_gene_name"), "ensembl_gene_id_version", transcript_ids, mart)  # create annotation to gene symbol key
  import_genes <- ann2gene[match(transcript_ids$V1, ann2gene$ensembl_gene_id_version),2]  # convert the list of annotations to list of corresponding gene symbols
  quant_array <- na.omit(t(read.csv(quant_csv, header = FALSE)))  # import the quantitation from alignment
  quant_array <- data.frame(t(quant_array))   # transpose dataframe to make more friendly to downstream manipulation
  colnames(quant_array) <- import_genes   # label quantitation array by gene name
  for (i in 1:nrow(quant_array)) {
    row.names(quant_array[i,]) <- paste0(batch_id,"_",sample_id,"_",i)
    print(paste0(batch_id,"_",sample_id,"_",i))
  }
  assign(paste0("quant_array_batch",batch_id,"_sample",sample_id), quant_array, env = .GlobalEnv)
}

# Create the array tallying metadata for each cell, only create after quant array has been made
AssembleMetadataArray <- function(batch_id, sample_id, region_id) {
  quant_reference <- get(paste0("quant_array_batch",batch_id,"_sample",sample_id), envir = .GlobalEnv)
  metadata_array <- data.frame(matrix(data = 0, nrow = nrow(quant_reference), ncol = 3))
  for (i in 1:nrow(quant_reference)) {
    row.names(metadata_array[i,]) <- paste0(batch_id,"_",sample_id,"_",i)
    print(paste0(batch_id,"_",sample_id,"_",i))
  }
  metadata_array[,1:3] <- c(batch_id, sample_id, region_id)
  colnames(metadata_array) <- c("Batch", "Sample", "Region")
  assign(paste0("metadata_array_batch",batch_id,"_sample",sample_id), metadata_array, env = .GlobalEnv)
}