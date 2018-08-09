# Matches cell IDs and annotations to values and assigns gene identities

source("setup.R")   # Import dependencies

# Create the quant array for gene expression by each cell and a corresponding array for metadata
AssembleArrays <- function(quant_csv, gene_ids, batch_id, sample_id, region_id, method_id) {
  transcript_ids <- read.table(gene_ids, colClasses = "character")  # import the annotations from alignment
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")   # connect to Ensembl via BioMart
  ann2gene <- getBM(c("ensembl_gene_id_version", "external_gene_name"), "ensembl_gene_id_version", transcript_ids, mart)  # create annotation to gene symbol key
  import_genes <- ann2gene[match(transcript_ids$V1, ann2gene$ensembl_gene_id_version),2]  # convert the list of annotations to list of corresponding gene symbols
  quant_array <- na.omit(t(read.csv(quant_csv, header = FALSE)))  # import the quantitation from alignment
  quant_array <- data.frame(t(quant_array))   # transpose dataframe to make more friendly to downstream manipulation
  colnames(quant_array) <- import_genes   # label quantitation array by gene name
  print("Quant array complete")
  metadata_array <- data.frame(matrix(data = 0, nrow = nrow(quant_array), ncol = 4))  # create separate array to store metadata for each cell
  colnames(metadata_array) <- c("batch", "sample", "region", "method")
  print("Metadata array complete")
  rownames(quant_array) <- paste0(batch_id,"_",sample_id,"_",seq(1:nrow(quant_array)))
  rownames(metadata_array) <- paste0(batch_id,"_",sample_id,"_",seq(1:nrow(metadata_array)))
  metadata_array$batch <- batch_id
  metadata_array$sample <- sample_id
  metadata_array$region <- region_id
  metadata_array$method <- method_id
  assign(paste0("quant_array_batch",batch_id,"_sample",sample_id), quant_array, env = .GlobalEnv)   # these need to be made like this so that it returns both with custom names
  assign(paste0("metadata_array_batch",batch_id,"_sample",sample_id), metadata_array, env = .GlobalEnv)
}