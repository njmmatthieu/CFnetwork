# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# # for get_network_nodes()
# source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/get_network_interactions.R")

# for endpoint_tag(), preprocess_PPI_network(), tag_weird_endpoints() and remove_weird_endpoints()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_visualization_helper.R")
# for remove_expression_interactions(), get_node_type()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/simplify_network_helper.R")

# for kegg_pathways_nodes.carac and effect_arrow.df
# source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_utilss.R")
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_nodes_carac.RData")

# To add an interactions
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_manual_curation.R")

# KEGG pathways from Omnipath
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_list.RData")
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/symbols_from_kegg_pathways_from_omnipathR_list.RData")

# for all_CFTR_interactors.PPI
source(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/extract_CyFi_MAP_CFTR_interactors.R")
# for extend_to_CFTR_interactors
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/CFTR_interactors_helper.R")

# Select pathway
pathway_name <- "Regulation of actin cytoskeleton"

RegActinCyto_kegg_pathway_df.interactions <- kegg_pathway_df.signaling.interactions.list.final[names(kegg_pathway_df.signaling.interactions.list.final)==pathway_name][[1]]

RegActinCyto_kegg_pathway_df.interactions.unique <- unique(RegActinCyto_kegg_pathway_df.interactions[,!(colnames(RegActinCyto_kegg_pathway_df.interactions) %in% c("source", 
                                                                                                                                                                   "target", 
                                                                                                                                                                   "relation_id", 
                                                                                                                                                                   "kegg_id_source",
                                                                                                                                                                   "kegg_id_target"))])

RegActinCyto_kegg_pathway_df.interactions.unique$manually_added <- FALSE

# Preprocess

# 
RegActinCyto_kegg_pathway_df.nodes <- kegg_pathway_df.signaling.nodes.list.final[names(kegg_pathway_df.signaling.nodes.list.final)==pathway_name][[1]]

# 
RegActinCyto_kegg_PPI_network <- new("PPI_network",
                                     interactions=RegActinCyto_kegg_pathway_df.interactions.unique,
                                     nodes=RegActinCyto_kegg_pathway_df.nodes)

## replace DDX58 by RIGI
RegActinCyto_kegg_PPI_network <- preprocess_PPI_network(RegActinCyto_kegg_PPI_network)

# RegActinCyto_kegg_PPI_network@nodes$pathway_name <- NULL
# Preprocess
kegg_pathways_corrections <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_diff_pathways_corrections_w_EZR_2022_11_28.txt",
                                        sep = "\t",
                                        header = T)
kegg_pathways_corrections.RegActinCyto <- kegg_pathways_corrections[which(kegg_pathways_corrections$pathway_name==pathway_name),]

RegActinCyto_kegg_PPI_network.corrected <- correct_PPI_network(PPI_network = RegActinCyto_kegg_PPI_network,
                                                               Symbol2UniprotID = kegg_pathways_nodes.carac.corrected,
                                                               corrections.df = kegg_pathways_corrections.RegActinCyto)

# Preprocess
RegActinCyto_kegg_PPI_network.corrected@interactions$interaction_id <- apply(X = RegActinCyto_kegg_PPI_network.corrected@interactions[,c("uniprot_source",
                                                                                                     "uniprot_target",
                                                                                                     "effect")],
                                                           MARGIN = 1,
                                                           FUN = function(x) {
                                                             # print(x)
                                                             return(paste(x, collapse = "_"))
                                                           })
RegActinCyto_kegg_PPI_network.corrected@interactions[,c("type", "pathway_name", "manually_added")] <- NULL


# Look if there are some CFTR interactors in the KEGG pathway
lapply(interactors.list, function(list_of_interactors){
  return(list_of_interactors[which(list_of_interactors %in% RegActinCyto_kegg_PPI_network.corrected@nodes$Symbol)])
})

# Extend CFTR interactions
RegActinCyto_kegg_PPI_network.CFTR_extended <- extend_to_CFTR_interactors(RegActinCyto_kegg_PPI_network.corrected)

# Node types 
# Endpoints: TF and Caspases (Apoptosis)
RegActinCyto_kegg_PPI_network.CFTR_extended <- endpoint_tag(RegActinCyto_kegg_PPI_network.CFTR_extended)
# Endpoints that are neither TF nor Capsapses
RegActinCyto_kegg_PPI_network.CFTR_extended <- tag_weird_endpoints(RegActinCyto_kegg_PPI_network.CFTR_extended)


# 1- Remove expression interactions (dorothea and kegg)
RegActinCyto_kegg_PPI_network.without_expression <- remove_expression_interactions(RegActinCyto_kegg_PPI_network.CFTR_extended)


# 3- Style - SHOULD PUT THE ONES FROM THE WHOLE NETWORK


RegActinCyto_kegg_PPI_network.without_expression@nodes <- merge(RegActinCyto_kegg_PPI_network.without_expression@nodes[,c("Symbol","pathway_name")],
                                                       CF_PPI_network.scc@nodes[,setdiff(colnames(CF_PPI_network.scc@nodes), diff_pathways)],
                                                       by = "Symbol",
                                                       all.x = T)

# RegActinCyto_kegg_PPI_network.tagged <- tag_prot_cat(RegActinCyto_kegg_PPI_network.without_expression)
# RegActinCyto_kegg_PPI_network.final <- get_node_type(RegActinCyto_kegg_PPI_network.tagged)

write.table(RegActinCyto_kegg_PPI_network.without_expression@interactions,
            file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                         pathway_name,
                         "_interactions_w_CFTR_interactors_2022_12_02.tsv",
                         sep = ""),
            sep = "\t",
            row.names = F,
            quote = FALSE)

write.table(RegActinCyto_kegg_PPI_network.without_expression@nodes,
            file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                         pathway_name,
                         "_nodes_w_CFTR_interactors_2022_12_02.tsv",
                         sep = ""),
            sep = "\t",
            row.names = F,
            quote = FALSE)
