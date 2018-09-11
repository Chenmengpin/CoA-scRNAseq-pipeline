# produces GO annotations for any group of genes we want to examine

# create list of GO terms
GetGOTerms <- function(gene_metadata) {
  gene_list <- rownames(gene_metadata)
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")   # connect to Ensembl via BioMart
  gene2GO <- getBM(attributes = c("external_gene_name", "name_1006"), filters = "external_gene_name", values = gene_list, mart = mart)  # create gene to GO term key
  GO_list_all <- dcast(data = gene2GO, formula = external_gene_name ~ name_1006, fun.aggregate = length)  # converts from melted to a binary array
  rownames(GO_list_all) <- GO_list_all$external_gene_name
  GO_list_all <- GO_list_all[, -(1:2)]  # this is needed to clean the array
  assign("GO_list", GO_list_all, env = .GlobalEnv)
}

WriteIDInferenceGenes <- function(GO_list) {
  IDInferenceGenes <- list()
  glutamate_genes <- rownames(GO_list[(GO_list$`adenylate cyclase inhibiting G-protein coupled glutamate receptor activity` | 
                                       GO_list$`adenylate cyclase-inhibiting G-protein coupled glutamate receptor signaling pathway` |
                                       GO_list$`AMPA glutamate receptor clustering` |
                                       GO_list$`asymmetric, glutamatergic, excitatory synapse` |
                                       GO_list$`G-protein coupled glutamate receptor signaling pathway` |
                                       GO_list$`glutamine biosynthetic process` |
                                       GO_list$`glutamate receptor activity` |
                                       GO_list$`glutamate secretion, neurotransmission` |
                                       GO_list$`glutamatergic postsynaptic density` |
                                       GO_list$`ionotropic glutamate receptor activity` |
                                       GO_list$`ionotropic glutamate receptor complex` |
                                       GO_list$`L-glutamate transmembrane transport` |
                                       GO_list$`NMDA glutamate receptor activity` |
                                       GO_list$`response to L-glutamate` |
                                       GO_list$`type 3 metabotropic glutamate receptor binding` == 1),])
  GABAglycine_genes <- rownames(GO_list[(GO_list$`cerebral cortex GABAergic interneuron fate commitment` | 
                                         GO_list$`G-protein coupled GABA receptor activity` |
                                         GO_list$`GABA-A receptor activity` |
                                         GO_list$`gamma-aminobutyric acid biosynthetic process` |
                                         GO_list$`gamma-aminobutyric acid receptor clustering` |
                                         GO_list$`gamma-aminobutyric acid transport` |
                                         GO_list$`glutamate decarboxylation to succinate` |
                                         GO_list$`extracellularly glycine-gated chloride channel activity` |
                                         GO_list$`glycine decarboxylation via glycine cleavage system` |
                                         GO_list$`glycine hydroxymethyltransferase activity` |
                                         GO_list$`glycine transmembrane transporter activity` |
                                         GO_list$`glycine transport` |
                                         GO_list$`synaptic transmission, glycinergic` == 1),])
  Neuropeptide_genes <- rownames(GO_list[(GO_list$`cholecystokinin receptor activity` | 
                                          GO_list$`neuropeptide binding` |
                                          GO_list$`neuropeptide hormone activity` |
                                          GO_list$`neuropeptide signaling pathway` |
                                          GO_list$`retrograde trans-synaptic signaling by neuropeptide, modulating synaptic transmission` |
                                          GO_list$`substance P catabolic process` |
                                          GO_list$`trans-synaptic signaling by trans-synaptic complex, modulating synaptic transmission` |
                                          GO_list$`tachykinin receptor signaling pathway` == 1),])
  Misc_genes <- rownames(GO_list[(GO_list$`cannabinoid receptor activity` | 
                                  GO_list$`regulation of endocannabinoid signaling pathway` |
                                  GO_list$`regulation of retrograde trans-synaptic signaling by endocanabinoid` |
                                  GO_list$`retrograde trans-synaptic signaling by endocannabinoid` |
                                  GO_list$`type 1 cannabinoid receptor binding` |
                                  GO_list$`histamine biosynthetic process` |
                                  GO_list$`histamine catabolic process` |
                                  GO_list$`histamine receptor activity` |
                                  GO_list$`histamine transport` |
                                  GO_list$`G-protein coupled adenosine receptor activity` |
                                  GO_list$`G-protein coupled purinergic nucleotide receptor activity` |
                                  GO_list$`purinergic nucleotide receptor activity` == 1),])
}