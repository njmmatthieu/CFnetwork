# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# for kegg_pathways_nodes.carac.corrected and effect_arrow.df
# source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_utils.R")

# for endpoint_tag(), 
# gene_symbol_sanity_check(), 
# tag_weird_endpoints(),
# remove_weird_endpoints(), 
# remove_expression_interactions(), 
# remove_indirect_interactions(), 
# remove_same_interactions(), 
# binding_interaction() 
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/simplify_network_helper.R")

# for get_node_type(),
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_visualization_helper.R")

# KEGG DIFF PATHWAYS - All proteins
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/kegg_diff_pathways_interactions_with_CFTR_interactors_df_2023_07_10.RData")
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/kegg_diff_pathways_nodes_with_CFTR_interactors_df_2023_07_10.RData")


# CF_PPI_network.CFTR_extended <- new("PPI_network",
#                                 interactions=kegg_diff_pathways_interactions,
#                                 nodes=kegg_diff_pathways_nodes)

CF_PPI_network.CFTR_extended <- new("PPI_network",
                                    interactions=CF_PPI_network.CFTR_extended.interactions,
                                    nodes=CF_PPI_network.CFTR_extended.nodes)

# ###
# OPTIONS - Extend CFTR interactions
# ###
# problem with kegg_pathways_nodes.carac.corrected (from kegg_diff_pathways.Rmd)
# problem with effect_arrow.df (from kegg_pathways_utils.R)
# CF_PPI_network.CFTR_extended <- extend_to_CFTR_interactors(CF_PPI_network)

# ### 
# A - OMNIPATH DB CURATIONS
# ###

# ### 
# 1 - REMOVE EXPRESSION INTERACTIONS
# ###

# Tag Endpoints: TF and Caspases (Apoptosis)
CF_PPI_network.CFTR_extended.endpoint_tag <- dorothea_tag(CF_PPI_network.CFTR_extended)


# Remove expression interactions (dorothea and kegg)
CF_PPI_network.without_expression <- remove_expression_interactions(CF_PPI_network.CFTR_extended.endpoint_tag)

# ### 
# 2 - REMOVE INDIRECT INTERACTIONS
# ###

CF_PPI_network.direct <- remove_indirect_interactions(CF_PPI_network.without_expression)

# # Tag receptors and receptor ligands
CF_PPI_network.direct.rep_tag <- tag_prot_cat(CF_PPI_network.direct)

CF_PPI_network.direct.node_type <- get_node_type(CF_PPI_network.direct.rep_tag,
                                              # interactors = CFTR_interactors,
                                              include_weird_endpoints = FALSE)

# ### 
# B - NETWORK CLEANING
# ###

# ### 
# 3 - REMOVE SAME INTERACTIONS: 1 phosphorylation (or ubiquitination) and 1 activation
# ###

CF_PPI_network.curated <- remove_same_interactions(CF_PPI_network.direct)



# ###
# 4 - REMOVE NON RECEPTORS THAT DON'T HAVE DOWNSTREAM INTERACTIONS
# ###

# Tag receptors and receptor ligands
CF_PPI_network.curated.rep_tag <- tag_prot_cat(CF_PPI_network.curated)

# Remove non receptor that don't have downstream interactions
# CF_PPI_network.lcc.curated.tagged <- tag_non_source_receptors_interactions(CF_PPI_network.lcc.curated)
CF_PPI_network.curated.2 <- remove_non_source_receptors(CF_PPI_network.curated.rep_tag)

CF_PPI_network.curated.2.node_type <- get_node_type(CF_PPI_network.curated.2,
                                                 # interactors = CFTR_interactors,
                                                 include_weird_endpoints = FALSE)

write.table(CF_PPI_network.curated.2.node_type@interactions,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df_2023_07_10.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)

write.table(CF_PPI_network.curated.2.node_type@nodes,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df_2022_07_10.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)


# ###
# 5 - SEEP SINGLE CONNECTED NETWORK
# ###

# CF_PPI_network.curated.tagged_lcc <- tag_largest_connected_component(CF_PPI_network.curated)
# pathway_components.count <- data.frame(table(CF_PPI_network.curated.tagged_lcc@nodes$pathway_components.membership))
# component_w_cascades <- pathway_components.count[which(pathway_components.count$Freq!=2),"Var1"]
# 
# for (component_id in component_w_cascades){
#   print(CF_PPI_network.curated.tagged_lcc@nodes[which(CF_PPI_network.curated.tagged_lcc@nodes$pathway_components.membership==component_id), "Symbol"])
# }
# CF_PPI_network.curated.tagged_lcc@nodes$pathway_components.membership <- as.factor(CF_PPI_network.curated.tagged_lcc@nodes$pathway_components.membership)
# 
# test <- CF_PPI_network.curated.tagged_lcc@nodes %>%
#   filter(pathway_components.membership %in% component_w_cascades) %>%
#   group_by(pathway_components.membership) %>% 
#   group_split() %>%
#   print

CF_PPI_network.lcc <- extract_largest_connected_component(CF_PPI_network.curated.2)
# CF_PPI_network.lcc.w_interactors <- extract_largest_connected_component(CF_PPI_network.curated.2.w_interactors)

# ###
# CHECKS
# ###

# ## Endpoints that are neither TF nor Capsapses
# CF_PPI_network.lcc.curated.2.endpoint_tag <- endpoint_tag(CF_PPI_network.lcc)
# 
# # Tag weird endpoints
# CF_PPI_network.lcc.curated.2.weird_endpoint <- tag_weird_endpoints(CF_PPI_network.lcc.curated.2.endpoint_tag)

# Network visualisation

# 1-  Node type

CF_PPI_network.lcc.node_type <- get_node_type(CF_PPI_network.lcc,
                                                # interactors = CFTR_interactors,
                                                include_weird_endpoints = FALSE)




write.table(CF_PPI_network.lcc.node_type@interactions,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_interactions_df_2023_07_10.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)

write.table(CF_PPI_network.lcc.node_type@nodes,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_nodes_df_2022_07_10.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)

# Analysis
# 
# library(gprofiler2)
# # KEGG signaling "gp__YJqt_TFVn_BE8"
# # KEGG signaling not expression gp__FZCS_noxZ_Snc
# upload_GMT_file(gmtfile = "/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_signaling_not_expression_from_omnipathR_gsea_2023_03_21.gmt")
# 
# CF_PPI_network.gost.res <- gost(query = CF_PPI_network.llc.curated.2.node_type@nodes$Symbol, 
#                                 organism = "gp__FZCS_noxZ_Snc", 
#                 ordered_query = FALSE,
#                 multi_query = FALSE, 
#                 significant = FALSE, 
#                 exclude_iea = FALSE, 
#                 measure_underrepresentation = FALSE, 
#                 evcodes = FALSE, 
#                 user_threshold = 0.05, 
#                 correction_method = "fdr", # change FDR
#                 # might change to all the genes in KEGG database
#                 # custom_bg = shared_genes_in_3_studies,
#                 # domain_scope = "custom", # change to custom to genes appearing in at least 3 studies
#                 numeric_ns = "", 
#                 sources = NULL, 
#                 as_short_link = FALSE)
# 
# CF_PPI_network.gost.res.df <- CF_PPI_network.gost.res$result
# CF_PPI_network.gost.res.df$ratio <- CF_PPI_network.gost.res.df$intersection_size/CF_PPI_network.gost.res.df$term_size   
# CF_PPI_network.gost.res.df$parents <- NULL

# write.table(CF_PPI_network.gost.res.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/KEGG_pathways_in_CF_network_2023_03_21.tsv",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)