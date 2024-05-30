# knitr::opts_knit$set(echo = TRUE, root.dir = normalizePath("../../"))

# for PPI_network class
source("scripts/pathways_to_network/network_utils.R")


# for endpoint_tag(), 
# gene_symbol_sanity_check(), 
# tag_weird_endpoints(),
# remove_weird_endpoints(), 
# remove_expression_interactions(), 
# remove_indirect_interactions(), 
# remove_same_interactions(), 
# binding_interaction() 
# tag_prot_cat()
source("scripts/network_analysis/simplify_network_helper.R")

# for get_node_type(),
source("scripts/network_analysis/network_visualization_helper.R")

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

# Import the data set in tidy tabular format
# NB: Multiple-value columns are kept as list-columns
(url <- latest_archive_url())
#> [1] "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
hgnc_dataset <- import_hgnc_dataset(url)
# Direct URL to the latest archive in TSV format

# # Reproducing Receptors
# receptors.df <- extract_receptors(CF_PPI_network.curated, 
#                                   genes_groups.df=hgnc_dataset)
  
PPI_network <- CF_PPI_network.curated
source_colname <- "genesymbol_source"
target_colname <- "genesymbol_target"
genes_groups.df <- hgnc_dataset


nodes_group <- genes_groups.df %>%
  dplyr::filter(symbol %in% PPI_network@nodes$Symbol) %>%
  dplyr::select(c('symbol', 
                  # 'status', 
                  'gene_group'))
test_receptors <- nodes_group %>%
  filter(symbol %in% c("PRLR", "EPOR", "GHR", "LEPR", "PLEKHG5", "CSF3R", "MPL"))
test_receptor_ligands <- nodes_group %>%
  filter(symbol %in% c("PRL", "EPO", "CSH1", "LEP", "TNFSF15", "CSF3"))

# Receptors
receptors.df <- nodes_group %>%
  filter(map_lgl(gene_group, ~  any(grepl(pattern = 'receptor', x = .x) &
                                      !(any(sapply(c('Receptor ligands',
                                                     'Estrogen receptors'),
                                                   function(fam) fam %in% .x))))))
  
  return(receptors.df)
  
}