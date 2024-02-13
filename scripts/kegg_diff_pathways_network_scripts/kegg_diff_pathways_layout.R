# knitr::opts_knit$set(echo = TRUE, root.dir = normalizePath("../../"))

# for PPI_network class
source("scripts/pathways_to_network/network_utils.R")

# for kegg_pathways_nodes.carac.corrected and effect_arrow.df
# source("scripts/kegg_pathways_scripts/kegg_pathways_utils.R")

# for endpoint_tag(), 
# gene_symbol_sanity_check(), 
# tag_weird_endpoints(),
# remove_weird_endpoints(), 
# remove_expression_interactions(), 
# remove_indirect_interactions(), 
# remove_same_interactions(), 
# binding_interaction() 
source("scripts/kegg_diff_pathways_network_scripts/simplify_network_helper.R")

# for get_node_type(),
source("scripts/kegg_diff_pathways_network_scripts/network_visualization_helper.R")

# KEGG DIFF PATHWAYS - All proteins
load("data/kegg_diff_pathways_network/kegg_diff_pathways_interactions_with_CFTR_interactors_df.RData")
load("data/kegg_diff_pathways_network/kegg_diff_pathways_nodes_with_CFTR_interactors_df.RData")

CF_PPI_network.CFTR_extended <- new("PPI_network",
                                    interactions=CF_PPI_network.CFTR_extended.interactions,
                                    nodes=CF_PPI_network.CFTR_extended.nodes)

# ### 
# A - OMNIPATH DB CURATIONS
# ###

# ### 
# 1 - REMOVE EXPRESSION INTERACTIONS
# ###

# Tag Endpoints: TF and Caspases (Apoptosis)
CF_PPI_network.CFTR_extended.endpoint_tag <- 
  dorothea_tag(CF_PPI_network.CFTR_extended)


# Remove expression interactions (dorothea and kegg)
CF_PPI_network.without_expression <- 
  remove_expression_interactions(CF_PPI_network.CFTR_extended.endpoint_tag)

# ### 
# 2 - REMOVE INDIRECT INTERACTIONS
# ###

CF_PPI_network.direct <- 
  remove_indirect_interactions(CF_PPI_network.without_expression)

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

# write.table(CF_PPI_network.curated.2.node_type@interactions,
#             file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(CF_PPI_network.curated.2.node_type@nodes,
#             file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)


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




# write.table(CF_PPI_network.lcc.node_type@interactions,
#             file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_interactions_df.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(CF_PPI_network.lcc.node_type@nodes,
#             file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_nodes_df.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)