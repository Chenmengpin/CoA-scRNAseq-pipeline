# produces GO annotations for any group of genes we want to examine

# create list of GO terms
GetGOTerms <- function(gene_metadata) {
  gene_list <- rownames(gene_metadata)
  print("Querying BioMart")
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")   # connect to Ensembl via BioMart
  gene2GO <- getBM(attributes = c("external_gene_name", "name_1006"), filters = "external_gene_name", values = gene_list, mart = mart)  # create gene to GO term key
  GO_library <- getBM(attributes = c("external_gene_name", "name_1006", "go_id", "definition_1006"), filters = "external_gene_name", values = gene_list, mart = mart)  # create gene to GO term key
  print("Building GO array")
  GO_array <- dcast(data = gene2GO, formula = external_gene_name ~ name_1006, fun.aggregate = length)  # converts from melted to a binary array
  rownames(GO_array) <- GO_array$external_gene_name
  print("Writing GO library and mappings")
  GO_library <- GO_library[GO_library$go_id != "",]   # removes empty elements before making mappings
  GO_mappings <- list()
  for (i in 1:length(gene_list)) {
    GO_mappings[[i]] <- GO_library$go_id[GO_library$external_gene_name %in% gene_list[i]]
    print(i)
  }
  names(GO_mappings) <- gene_list
  GO_array <- GO_array[, -(1:2)]  # these are needed to clean the array
  GO_library <- GO_library[, -1]
  GO_library <- GO_library[!duplicated(GO_library$name_1006),]   # needed to remove duplicates from the master library
  colnames(GO_library) <- c("GO_Term", "GO_ID", "GO_Description")
  Glutamate_terms <- c("adenylate cyclase inhibiting G-protein coupled glutamate receptor activity", "adenylate cyclase-inhibiting G-protein coupled glutamate receptor signaling pathway", "AMPA glutamate receptor clustering",
                       "asymmetric, glutamatergic, excitatory synapse", "G-protein coupled glutamate receptor signaling pathway", "glutamine biosynthetic process", "glutamate receptor activity",
                       "glutamate secretion, neurotransmission", "glutamatergic postsynaptic density", "ionotropic glutamate receptor activity", "ionotropic glutamate receptor complex", "L-glutamate transmembrane transport",
                       "NMDA glutamate receptor activity", "response to L-glutamate", "type 3 metabotropic glutamate receptor binding")
  GABAglycine_terms <- c("cerebral cortex GABAergic interneuron fate commitment", "G-protein coupled GABA receptor activity", "GABA-A receptor activity", "gamma-aminobutyric acid biosynthetic process", "glycine transport",
                         "gamma-aminobutyric acid receptor clustering", "gamma-aminobutyric acid transport", "glutamate decarboxylation to succinate", "extracellularly glycine-gated chloride channel activity",
                         "glycine decarboxylation via glycine cleavage system", "glycine hydroxymethyltransferase activity", "glycine transmembrane transporter activity", "synaptic transmission, glycinergic")
  Monoamine_terms <- c("monoamine transport", "regulation of serotonin uptake", "isoquinoline alkaloid metabolic process", "dopamine catabolic process", "dopamine metabolic process", "dopamine neurotransmitter receptor activity",
                       "dopamine uptake involved in synaptic transmission", "positive regulation of dopamine metabolic process", "adrenergic receptor activity", "alpha-1A adrenergic receptor binding",
                       "alpha-2A adrenergic receptor binding", "beta-3 adrenergic receptor binding", "noradrenergic neuron development", "histamine biosynthetic process", "histamine catabolic process", "histamine transport",
                       "serotonin binding", "serotonin biosynthetic process", "serotonin transmembrane transporter activity", "negative regulation of serotonin secretion")
  MiscNT_terms <- c("cannabinoid receptor activity", "regulation of endocannabinoid signaling pathway", "regulation of retrograde trans-synaptic signaling by endocanabinoid", "retrograde trans-synaptic signaling by endocannabinoid",
                    "type 1 cannabinoid receptor binding", "acetylcholine binding", "acetylcholine receptor inhibitor activity", "acetylcholine transmembrane transporter activity", 
                    "acetylcholine-gated cation-selective channel activity", "acetylcholine-gated channel complex", "acetylcholinesterase activity", "choline transmembrane transporter activity", 
                    "negative regulation of skeletal muscle acetylcholine-gated channel clustering", "skeletal muscle acetylcholine-gated channel clustering", "synaptic transmission, cholinergic", 
                    "G-protein coupled adenosine receptor activity", "G-protein coupled purinergic nucleotide receptor activity", "purinergic nucleotide receptor activity")
  IonChannel_terms <- c("A-type (transient outward) potassium channel activity", "ATP-activated inward rectifier potassium channel activity", "inward rectifier potassium channel activity",
                        "L-type voltage-gated calcium channel complex", "open rectifier potassium channel activity", "outward rectifier potassium channel activity", "voltage-gated anion channel activity",
                        "voltage-gated calcium channel complex", "voltage-gated chloride channel activity", "voltage-gated ion channel activity", "voltage-gated potassium channel activity", "voltage-gated sodium channel complex")
  PeptideHormoneNT_terms <- c("cholecystokinin receptor activity", "neuropeptide binding", "neuropeptide hormone activity", "retrograde trans-synaptic signaling by neuropeptide, modulating synaptic transmission", 
                              "neuropeptide signaling pathway", "substance P catabolic process", "trans-synaptic signaling by trans-synaptic complex, modulating synaptic transmission", "tachykinin receptor signaling pathway",
                              "hormone activity", "hormone binding", "hormone receptor binding", "hormone-mediated signaling pathway", "peptide hormone binding", "peptide hormone receptor binding", "protein-hormone receptor activity")
  Adhesion_terms <- c("cadherin binding involved in cell-cell adhesion", "calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules", "calcium-independent cell-cell adhesion via plasma membrane cell-adhesion molecules",
                      "calcium-independent cell-matrix adhesion", "cell adhesion mediated by integrin", "cell adhesion molecule binding", "cell-cell adhesion mediated by cadherin", "cell-cell adhesion mediated by integrin",
                      "cell-cell adhesion involved in neuronal-glial interactions involved in cerebral cortex radial glia guided migration", "cell-cell adhesion mediator activity", "cell-cell adhesion via plasma-membrane adhesion molecules",
                      "collagen binding involved in cell-matrix adhesion", "desmosome", "heterophilic cell-cell adhesion via plasma membrane cell adhesion molecules", "heterotypic cell-cell adhesion", "homotypic cell-cell adhesion",
                      "homophilic cell adhesion via plasma membrane adhesion molecules", "integrin binding involved in cell-matrix adhesion", "integrin complex", "neuron cell-cell adhesion", "synaptic membrane adhesion",
                      "neuronal-glial interaction involved in cerebral cortex radial glia guided migration")
  Pathfinding_terms <- c("anterior/posterior axon guidance", "axon extension involved in axon guidance", "axon guidance receptor activity", "commissural neuron axon guidance", "dorsal spinal cord interneuron anterior axon guidance",
                         "motor neuron axon guidance", "negative regulation of axon extension involved in axon guidance", "netrin-activated signaling pathway", "olfactory bulb axon guidance", "neuron projection guidance",
                         "planar cell polarity pathway involved in axon guidance", "positive regulation of axon extension involved in axon guidance", "positive regulation of retinal ganglion cell axon guidance", "neuropilin binding",
                         "regulation of axon guidance", "retinal ganglion cell axon guidance", "sensory neuron axon guidance", "semaphorin-plexin signaling pathway", "semaphorin receptor binding", "neuron projection maintenance")
  Chemosensory_terms <- c("olfactory receptor activity", "pheromone receptor activity", "trace-amine receptor activity", "bitter taste receptor activity", "sensory perception of sour taste", "sensory perception of umami taste",
                          "taste receptor activity")
  print("Sorting genes into functional groups")
  GO_Term_List <- list(Glutamate_terms, GABAglycine_terms, Monoamine_terms, MiscNT_terms, IonChannel_terms, PeptideHormoneNT_terms, Adhesion_terms, Pathfinding_terms, Chemosensory_terms)  # group all the known terms together
  GO_IDs <- list()  # needed to prevent the loop from crashing
  all_genes <- c("Initial_entry")
  for (i in 1:9){
    GO_annotation <- GO_library[GO_library$GO_Term %in% GO_Term_List[[i]],]
    geneIDs_GO <- which(apply(GO_array[, colnames(GO_array) %in% GO_Term_List[[i]]] == 1, 1, any))  # find the row IDs of the chosen genes
    genes_GO <- rownames(GO_array[geneIDs_GO, ])
    genes_GO <- genes_GO[!all_genes %in% genes_GO]
    GO_IDs[[i]] <- list(genes_GO, GO_annotation)
    names(GO_IDs[[i]]) <- c("Genes", "GO_Annotation")
    all_genes <- c(all_genes, genes_GO)
    print(i)
  }
  names(GO_IDs) <- c("Glutamate", "GABA_glycine", "Monoamine", "Misc_NT", "IonChannel", "Peptide_Hormone_NT", "Adhesion", "Pathfinding", "Chemosensation")
  assign("GO_array", GO_array, env = .GlobalEnv)
  assign("GO_library", GO_library, env = .GlobalEnv)
  assign("GO_IDs", GO_IDs, env = .GlobalEnv)
  assign("GO_mappings", GO_mappings, env = .GlobalEnv)
}

GetEnrichedGO <- function(GO_genes, GO_mappings, GO_name) {
  Enriched_GOs <- list()
  for (i in 1:length(GO_genes)) {
    GO_data_BP <- new("topGOdata", ontology = "BP", allGenes = GO_genes[[i]],
                      annot = annFUN.gene2GO, gene2GO = GO_mappings)
    GO_data_MF <- new("topGOdata", ontology = "MF", allGenes = GO_genes[[i]],
                      annot = annFUN.gene2GO, gene2GO = GO_mappings)
    GO_data_CC <- new("topGOdata", ontology = "CC", allGenes = GO_genes[[i]],
                      annot = annFUN.gene2GO, gene2GO = GO_mappings)
    GOtest <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher", allMembers = names(GO_genes[[i]]),
                  sigMembers = names(GO_genes[[i]][GO_genes[[i]] == 1]))
    resultFisher_BP <- getSigGroups(GO_data_BP, test.stat = GOtest)
    resultTable_BP <- GenTable(GO_data_BP, classic = resultFisher_BP, orderBy = "classic", 
                            ranksOf = "classic", topNodes = length(usedGO(GO_data_BP)))
    resultTable_BP$classic <- p.adjust(resultTable_BP$classic, method = "BH")
    resultFisher_MF <- getSigGroups(GO_data_MF, test.stat = GOtest)
    resultTable_MF <- GenTable(GO_data_MF, classic = resultFisher_MF, orderBy = "classic", 
                               ranksOf = "classic", topNodes = length(usedGO(GO_data_MF)))
    resultTable_MF$classic <- p.adjust(resultTable_MF$classic, method = "BH")
    resultFisher_CC <- getSigGroups(GO_data_CC, test.stat = GOtest)
    resultTable_CC <- GenTable(GO_data_CC, classic = resultFisher_CC, orderBy = "classic", 
                               ranksOf = "classic", topNodes = length(usedGO(GO_data_CC)))
    resultTable_CC$classic <- p.adjust(resultTable_CC$classic, method = "BH")
    enrichedGO <- rbind.data.frame(resultTable_BP[resultTable_BP$classic < .05,],
                                   resultTable_MF[resultTable_MF$classic < .05,],
                                   resultTable_CC[resultTable_CC$classic < .05,])
    Enriched_GOs[[i]] <- enrichedGO
  }
  names(Enriched_GOs) <- names(GO_genes)
  assign(GO_name, Enriched_GOs, env = .GlobalEnv)
}
