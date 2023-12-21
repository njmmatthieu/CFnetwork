# conda env : doro-r4

library(OmnipathR)
library(data.table)
library(tidyverse)
library(dorothea)

kegg_pathways <- data.frame(kegg_pathway_list())

kegg_pathways.signalling <- kegg_pathways[which(grepl("hsa04[[:digit:]]{3}",kegg_pathways$id)),]

# Remove pathways describing diseases from the pathways list
disease_signaling.ids <- c("hsa04930",
                           "hsa04931",
                           "hsa04932",
                           "hsa04933",
                           "hsa04934",
                           "hsa04936",
                           "hsa04940",
                           "hsa04950")
kegg_pathways.signalling <- kegg_pathways.signalling[which(!kegg_pathways.signalling$id %in% disease_signaling.ids),]

# GMT

## All genes

kegg_pathways.signaling.list <- lapply(1:dim(kegg_pathways.signalling)[1], 
                                       function(pathway_i){
  
  # to test
  # pathway_i <- 5
  kegg_pathways.signalling$name[pathway_i]
  
  kegg_pathway_df <- kegg_pathway_download(kegg_pathways.signalling$id[pathway_i], 
                                           process = TRUE)
  kegg_pathway_genes <- unique(c(kegg_pathway_df$genesymbol_source,
                                 kegg_pathway_df$genesymbol_target))
  # print(kegg_pathway_genes)
  return(as.data.frame(t(as.data.frame(c(kegg_pathways.signalling$name[pathway_i], 
           'KEGG',
           kegg_pathway_genes), ))))
})


kegg_pathways.signaling.df <- rbindlist(kegg_pathways.signaling.list,
                                        fill = T)

# write.table(kegg_pathways.signaling.df,
#             file = "kegg_pathways/kegg_from_omnipathR_gsea_2022_09_07.gmt",
#             sep = "\t",
#             row.names = FALSE,
#             col.names = FALSE,
#             na = "",
#             quote = FALSE)

# Interactions data frame

kegg_pathway_df.signaling.interactions.list <- lapply(1:dim(kegg_pathways.signalling)[1], 
                                                      function(pathway_i){
  
  pathway_name <- kegg_pathways.signalling$name[pathway_i]
  
  kegg_pathways.signalling.df <- data.frame(kegg_pathway_download(kegg_pathways.signalling$id[pathway_i], 
                                                                  process = TRUE))
  kegg_pathways.signalling.df$pathway_name <- rep(pathway_name,
                                                  dim(kegg_pathways.signalling.df)[1])
  
  return(kegg_pathways.signalling.df)
  
})
names(kegg_pathway_df.signaling.interactions.list) <- lapply(kegg_pathway_df.signaling.interactions.list,function(df){
  return(unique(df$pathway_name))
})
kegg_pathway_df.signaling.interactions.list.final <- 
  kegg_pathway_df.signaling.interactions.list[lapply(kegg_pathway_df.signaling.interactions.list, function(df) {
    return(dim(df)[2])})==13]

# save(kegg_pathway_df.signaling.interactions.list.final,
#      file = "kegg_pathways/kegg_pathways_from_omnipath_list.RData")

kegg_pathway_df.signaling.interactions.df <- rbindlist(kegg_pathway_df.signaling.interactions.list.final)

# write.table(kegg_pathway_df.signaling.interactions.df,
#             file = "kegg_pathways/kegg_pathways_from_omnipathR.txt",
#             sep = "\t",
#             # row.names = FALSE,
#             # col.names = FALSE,
#             # na = "",
#             quote = FALSE)

# Nodes data frame

kegg_pathway_df.signaling.nodes.list <- lapply(1:dim(kegg_pathways.signalling)[1], function(pathway_i){
  
  pathway_name <- kegg_pathways.signalling$name[pathway_i]
  
  kegg_pathway_df <- data.frame(kegg_pathway_download(kegg_pathways.signalling$id[pathway_i], 
                                                      process = TRUE))
  kegg_pathway_symbols_df <- data.frame(Symbol = unique(c(kegg_pathway_df$genesymbol_source,
                                                           kegg_pathway_df$genesymbol_target)))
  kegg_pathway_symbols_df$pathway_name <- rep(pathway_name,
                                              dim(kegg_pathway_symbols_df)[1])
  
  return(kegg_pathway_symbols_df)
  
})
names(kegg_pathway_df.signaling.nodes.list) <- lapply(kegg_pathway_df.signaling.nodes.list,function(df){
  return(unique(df$pathway_name))
})

kegg_pathway_df.signaling.nodes.list.final <- 
  kegg_pathway_df.signaling.nodes.list[lapply(kegg_pathway_df.signaling.nodes.list, function(df) {
    return(dim(df)[2])})==2]

# save(kegg_pathway_df.signaling.nodes.list.final,
#      file = "kegg_pathways/symbols_from_kegg_pathways_from_omnipathR_list.RData")

kegg_pathway_df.signaling.nodes.df <- rbindlist(kegg_pathway_df.signaling.nodes.list.final)

# write.table(kegg_pathway_df.signaling.nodes.df,
#             file = "kegg_pathways/symbols_from_kegg_pathways_from_omnipathR.txt",
#             sep = "\t",
#             # row.names = FALSE,
#             # col.names = FALSE,
#             # na = "",
#             quote = FALSE)
