# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# for kegg_pathways_nodes.carac.corrected and effect_arrow.df
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_utils.R")

# To add an interactions
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_manual_curation.R")

# for endpoint_tag(), 
# gene_symbol_sanity_check(), 
# tag_weird_endpoints(),
# remove_weird_endpoints(), 
# remove_expression_interactions(), 
# remove_indirect_interactions(), 
# remove_same_interactions(), 
# binding_interaction() 
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/simplify_network_helper.R")


# for all_CFTR_interactors.PPI
source(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/extract_CyFi_MAP_CFTR_interactors.R")
# extend to CFTR interactors
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/CFTR_interactors_helper.R")

# for get_node_type(),
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_visualization_helper.R")


# Select pathway
pathway_name <- "TNF signaling pathway"

TNF_kegg_pathway_df.interactions <- kegg_pathway_df.signaling.interactions.list.final[names(kegg_pathway_df.signaling.interactions.list.final)==pathway_name][[1]]

TNF_kegg_pathway_df.interactions.unique <- unique(TNF_kegg_pathway_df.interactions[,!(colnames(TNF_kegg_pathway_df.interactions) %in% c("source", 
                                                                                                                                                                   "target", 
                                                                                                                                                                   "relation_id", 
                                                                                                                                                                   "kegg_id_source",
                                                                                                                                                                   "kegg_id_target"))])

TNF_kegg_pathway_df.interactions.unique$manually_added <- FALSE

# Preprocess

# 
TNF_kegg_pathway_df.nodes <- kegg_pathway_df.signaling.nodes.list.final[names(kegg_pathway_df.signaling.nodes.list.final)==pathway_name][[1]]

# 
TNF_kegg_PPI_network <- new("PPI_network",
                                     interactions=TNF_kegg_pathway_df.interactions.unique,
                                     nodes=TNF_kegg_pathway_df.nodes)

## replace DDX58 by RIGI
TNF_kegg_PPI_network <- preprocess_PPI_network(TNF_kegg_PPI_network)

# TNF_kegg_PPI_network@nodes$pathway_name <- NULL
# Preprocess
kegg_pathways_corrections <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_diff_pathways_corrections_w_EZR_2022_11_28.txt",
                                        sep = "\t",
                                        header = T)
kegg_pathways_corrections.TNF <- kegg_pathways_corrections[which(kegg_pathways_corrections$pathway_name==pathway_name),]

TNF_kegg_PPI_network.corrected <- correct_PPI_network(PPI_network = TNF_kegg_PPI_network,
                                                               Symbol2UniprotID = kegg_pathways_nodes.carac.corrected,
                                                               corrections.df = kegg_pathways_corrections.TNF)

# Preprocess
TNF_kegg_PPI_network.corrected@interactions$interaction_id <- apply(X = TNF_kegg_PPI_network.corrected@interactions[,c("uniprot_source",
                                                                                                                                         "uniprot_target",
                                                                                                                                         "effect")],
                                                                             MARGIN = 1,
                                                                             FUN = function(x) {
                                                                               # print(x)
                                                                               return(paste(x, collapse = "_"))
                                                                             })
TNF_kegg_PPI_network.corrected@interactions[,c("type", "pathway_name", "manually_added")] <- NULL

# Look if there are some CFTR interactors in the KEGG pathway
lapply(interactors.list, function(list_of_interactors){
  return(list_of_interactors[which(list_of_interactors %in% TNF_kegg_PPI_network.corrected@nodes$Symbol)])
})

# Extend CFTR interactions
TNF_kegg_PPI_network.CFTR_extended <- extend_to_CFTR_interactors(TNF_kegg_PPI_network.corrected,
                                                                 include_CFTR = TRUE)
# No interactor to connect to this pathway

# Node types 
# Endpoints: TF and Caspases (Apoptosis)
TNF_kegg_PPI_network.CFTR_extended <- endpoint_tag(TNF_kegg_PPI_network.corrected)
# Endpoints that are neither TF nor Capsapses
TNF_kegg_PPI_network.CFTR_extended <- tag_weird_endpoints(TNF_kegg_PPI_network.CFTR_extended)


# 1- Remove expression interactions (dorothea and kegg)
TNF_kegg_PPI_network.without_expression <- remove_expression_interactions(TNF_kegg_PPI_network.CFTR_extended)

# 2- Remove indirect interactions
TNF_kegg_PPI_network.final <- remove_indirect_interactions(TNF_kegg_PPI_network.without_expression)

# 3- Style - SHOULD PUT THE ONES FROM THE WHOLE NETWORK
TNF_kegg_PPI_network.final@nodes <- merge(TNF_kegg_PPI_network.final@nodes[,c("Symbol","pathway_name")],
                                                                CF_PPI_network.scc@nodes[,setdiff(colnames(CF_PPI_network.scc@nodes), diff_pathways)],
                                                                by = "Symbol",
                                                                all.x = T)

# TNF_kegg_PPI_network.tagged <- tag_prot_cat(TNF_kegg_PPI_network.without_expression)
# TNF_kegg_PPI_network.final <- get_node_type(TNF_kegg_PPI_network.tagged)

write.table(TNF_kegg_PPI_network.final@interactions,
            file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                         pathway_name,
                         "_interactions_w_CFTR_interactors_2022_12_09.tsv",
                         sep = ""),
            sep = "\t",
            row.names = F,
            quote = FALSE)

write.table(TNF_kegg_PPI_network.final@nodes,
            file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                         pathway_name,
                         "_nodes_w_CFTR_interactors_2022_12_09.tsv",
                         sep = ""),
            sep = "\t",
            row.names = F,
            quote = FALSE)
