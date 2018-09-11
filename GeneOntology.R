# produces GO annotations for any group of genes we want to examine

# create list of GO terms
GetGOTerms <- function(gene_metadata) {
  gene_list <- rownames(gene_metadata)
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")   # connect to Ensembl via BioMart
  gene2GO <- getBM(attributes = c("external_gene_name", "name_1006"), filters = "external_gene_name", values = gene_list, mart = mart)  # create annotation to gene symbol key
  GO_list_all <- dcast(data = gene2GO, formula = external_gene_name ~ name_1006, fun.aggregate = length)
  assign("GO_list_all", GO_list_all, env = .GlobalEnv)
}