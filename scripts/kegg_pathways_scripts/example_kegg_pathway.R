# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# KEGG pathways from Omnipath
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_list.RData")
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_nodes_carac.RData")
kegg_pathway_df.signaling.interactions.df <- do.call("rbind", kegg_pathway_df.signaling.interactions.list.final)
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/symbols_from_kegg_pathways_from_omnipathR_list.RData")
kegg_pathway_df.signaling.nodes.df <- do.call("rbind", kegg_pathway_df.signaling.nodes.list.final)

# for corrections 
kegg_pathways_corrections <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_diff_pathways_corrections_w_EZR_2023_06_07.txt",
                                        sep = "\t",
                                        header = T)

# # To add an interactions
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_manual_curation.R")
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_utils.R")

# for remove_expression_interactions(), get_node_type()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/simplify_network_helper.R")

# for endpoint_tag(), preprocess_PPI_network(), tag_weird_endpoints() and remove_weird_endpoints()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_visualization_helper.R")

# # for extend_to_CFTR_interactors
# source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/CFTR_interactors_helper.R")


# KEGG pathways from Omnipath
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_list.RData")
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/symbols_from_kegg_pathways_from_omnipathR_list.RData")

# Select pathway
# example_pathway_name <- "NOD-like receptor signaling pathway"
# example_pathway_name <- "Viral protein interaction with cytokine and cytokine receptor"
# example_pathway_name <- "C-type lectin receptor signaling pathway"
# example_pathway_name <- "Cytosolic DNA-sensing pathway"
# example_pathway_name <- "MAPK signaling pathway" 
# example_pathway_name <- "RIG-I-like receptor signaling pathway"
# example_pathway_name <- "TNF signaling pathway"
# example_pathway_name <- "NF-kappa B signaling pathway"
# example_pathway_name <- "Toll-like receptor signaling pathway"
# example_pathway_name <- "IL-17 signaling pathway"
# example_pathway_name <- "Th17 cell differentiation"
# example_pathway_name <- "Cytokine-cytokine receptor interaction"
# example_pathway_name <- "T cell receptor signaling pathway"
example_pathway_name <- "Estrogen signaling pathway"
# example_pathway_name <- "Osteoclast differentiation"
# example_pathway_name <- "Regulation of actin cytoskeleton"

example_kegg_pathway_df.interactions <- kegg_pathway_df.signaling.interactions.list.final[names(kegg_pathway_df.signaling.interactions.list.final)==example_pathway_name][[1]]
# example_kegg_pathway_df.interactions$interaction_id <- apply(X = example_kegg_pathway_df.interactions[,c("uniprot_source",
#                                                                                      "uniprot_target",
#                                                                                      "effect")],
#                                                    MARGIN = 1,
#                                                    FUN = function(x) {
#                                                      # print(x)
#                                                      return(paste(x, collapse = "_"))
#                                                    })
example_kegg_pathway_df.interactions.unique <- unique(example_kegg_pathway_df.interactions[,!(colnames(example_kegg_pathway_df.interactions) %in% c("source", 
                                                                                                                      "target", 
                                                                                                                      "relation_id", 
                                                                                                                      "kegg_id_source",
                                                                                                                      "kegg_id_target"))])

example_kegg_pathway_df.interactions.unique$manually_added <- FALSE

# 
example_kegg_pathway_df.nodes <- kegg_pathway_df.signaling.nodes.list.final[names(kegg_pathway_df.signaling.nodes.list.final)==example_pathway_name][[1]]
# 
example_kegg_PPI_network <- new("PPI_network",
                                interactions=example_kegg_pathway_df.interactions.unique,
                                nodes=example_kegg_pathway_df.nodes)

## replace DDX58 by RIGI
example_kegg_PPI_network <- preprocess_PPI_network(example_kegg_PPI_network)


# Corrections preprocess
# kegg_pathways_corrections <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_diff_pathways_corrections_w_EZR_2022_12_22.txt",
#                                         sep = "\t",
#                                         header = T)
kegg_pathways_corrections.example <- kegg_pathways_corrections[which(kegg_pathways_corrections$pathway_name==example_pathway_name),]

example_kegg_PPI_network.corrected <- correct_PPI_network(PPI_network = example_kegg_PPI_network,
                                                                Symbol2UniprotID = kegg_pathways_nodes.carac.corrected,
                                                                corrections.df = kegg_pathways_corrections.example)

# Preprocess
example_kegg_PPI_network.corrected@interactions$interaction_id <- apply(X = example_kegg_PPI_network.corrected@interactions[,c("uniprot_source",
                                                                                                                                         "uniprot_target",
                                                                                                                                         "effect")],
                                                                             MARGIN = 1,
                                                                             FUN = function(x) {
                                                                               # print(x)
                                                                               return(paste(x, collapse = "_"))
                                                                             })
example_kegg_PPI_network.corrected@interactions[,c("type", "pathway_name", "manually_added")] <- NULL

# Extend CFTR interactions
# example_kegg_PPI_network.CFTR_extended <- extend_to_CFTR_interactors(example_kegg_PPI_network.corrected)

# Tag Endpoints: TF and Caspases (Apoptosis)
example_kegg_PPI_network.endpoint_tag <- dorothea_tag(example_kegg_PPI_network.corrected)


# # Remove expression interactions (dorothea and kegg)
# example_kegg_PPI_network.without_expression <- remove_expression_interactions(example_kegg_PPI_network.CFTR_extended.endpoint_tag)

# Tag receptors and receptor ligands
example_kegg_PPI_network.rep_tag <- tag_prot_cat(example_kegg_PPI_network.endpoint_tag)

# Node type
example_kegg_PPI_network.node_type <- get_node_type(example_kegg_PPI_network.rep_tag,
                                              # interactors = CFTR_interactors,
                                              include_weird_endpoints = FALSE)

# for diff_pathways
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/Transcriptomic studies/fgsea_comparison/fgsea_nes_diff_pathways_2022_09_15.RData")
diff_pathways <- rownames(fgsea_es_diff_pathways)
diff_pathways <- diff_pathways[which(diff_pathways!="AGE-RAGE signaling pathway in diabetic complications")]

# Loading leading.edges columns
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/kegg_diff_pathways_nodes_df_2023_04_06.RData")
leading.edges.columns <- colnames(kegg_diff_pathways_nodes)[grepl("leadingEdge", colnames(kegg_diff_pathways_nodes))]
leading.edges.column.pathway <- leading.edges.columns[grepl(example_pathway_name, leading.edges.columns)]

example_kegg_PPI_network.node_type@nodes <- merge(example_kegg_PPI_network.node_type@nodes,
                                        kegg_diff_pathways_nodes[,c("Symbol", leading.edges.column.pathway)],
                                        by = "Symbol",
                                        all.x = T)

write.table(example_kegg_PPI_network.node_type@interactions,
             file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                          example_pathway_name,
                          "_interactions_corrected_with_expression_2023_06_06.tsv",
                          sep = ""),
             sep = "\t",
             row.names = F,
             quote = FALSE)

write.table(example_kegg_PPI_network.node_type@nodes,
            file = paste("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/",
                         example_pathway_name,
                         "_nodes_corrected_with_expression_2023_06_06.tsv",
                         sep = ""),
            sep = "\t",
            row.names = F,
            quote = FALSE)

