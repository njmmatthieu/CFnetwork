# conda env : deseq

# library(hgnc)

# # Date of HGNC last update
# last_update()
# #> [1] "2022-04-04 02:37:12 UTC"
# 
# # Direct URL to the latest archive in TSV format
# (url <- latest_archive_url())
# #> [1] "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
# 
# # Import the data set in tidy tabular format
# # NB: Multiple-value columns are kept as list-columns
# hgnc_dataset <- import_hgnc_dataset(url)

# # KEGG pathways from Omnipath
# load("kegg_pathways_from_omnipath_list.RData")
# kegg_pathway_df.signaling.interactions.df <- do.call("rbind", kegg_pathway_df.signaling.interactions.list.final)
# load("symbols_from_kegg_pathways_from_omnipathR_list.RData")
# kegg_pathway_df.signaling.nodes.df <- do.call("rbind", kegg_pathway_df.signaling.nodes.list.final)

# kegg_pathways_nodes.carac <- hgnc_dataset %>%
#   dplyr::filter(symbol %in% unique(kegg_pathway_df.signaling.nodes.df$Symbol)) %>%
#   dplyr::select(c('symbol', 
#                   'uniprot_ids',
#                   'gene_group'))
# 
# kegg_pathways_nodes.carac[which(kegg_pathways_nodes.carac$symbol=="GNAS"),"uniprot_ids"][[1]] <- list("P63092")
# 
# save(kegg_pathways_nodes.carac,
#      file = "kegg_pathways_from_omnipath_nodes_carac.RData")

# Manual curation: for add_new_nodes_uniprot_ids
# source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_pathways_manual_curation.R")
# kegg_pathways_corrections <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/kegg_diff_pathways_corrections_w_EZR_2022_12_22.txt",
#                                         sep = "\t",
#                                         header = T)
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_pathways_from_omnipath_nodes_carac.RData")
kegg_pathways_nodes.carac.corrected <- add_new_nodes_uniprot_ids(Symbol2UniprotID = kegg_pathways_nodes.carac,
                                                                 corrections.df = kegg_pathways_corrections)


effect_arrow.df <- unique(kegg_pathway_df.signaling.interactions.df[,c("effect", "arrow")])
effect_arrow.df[which(effect_arrow.df$effect=="compound"),"effect"] <- paste(
  effect_arrow.df[which(effect_arrow.df$effect=="compound"),"effect"],
  effect_arrow.df[which(effect_arrow.df$effect=="compound"),"arrow"],
  sep = " - "
)

compound.df <- data.frame(effect=c(37, 89),
                          compound=c("PIP2", "IP3"))
