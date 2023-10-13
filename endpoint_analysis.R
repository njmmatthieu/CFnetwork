library(dorothea)
library(hgnc)
library(pheatmap)
library(tidyverse)

# for PPI_network class, get_network_nodes()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

####

# # to test
source_colname="genesymbol_source"
target_colname="genesymbol_target"
PPI_network <- CF_PPI_network.lcc.node_type

# # Endpoints: TF and Caspases (Apoptosis)
endpoint_tag <-function(PPI_network,
                        source_colname="genesymbol_source",
                        target_colname="genesymbol_target"){
  
  # # Caspases
  # apoptosis_nodes <- PPI_network@nodes$Symbol[grep("CASP", PPI_network@nodes$Symbol)]
  # PPI_network@nodes$endpoint_to_keep <- PPI_network@nodes$dorothea_tf | 
  #   PPI_network@nodes$Symbol %in% apoptosis_nodes
  # endpoint_to_keep.symbols <- PPI_network@nodes$Symbol[which(PPI_network@nodes$endpoint_to_keep)]
  #   # endpoints with binding interactions
  # endpoints_to_keep.non_binding <- sapply(endpoint_to_keep.symbols$Symbol, function(symbol){
  #   
  #   # to test
  #   symbol <- "STAT1"
  #   interactions.df <- PPI_network@interactions[which(PPI_network@interactions$genesymbol_source ==symbol |
  #                                                       PPI_network@interactions$genesymbol_target==symbol),]
  #   # remove binding interactions
  #   interactions.df <- interactions.df[which(interactions.df$arrow!="---"),]
  #   
  #   endpoints.targets <- setdiff(interactions.df[,target_colname], 
  #                                interactions.df[,source_colname])
  #   endpoints.targets <- setdiff(endpoints.targets, c(symbol))
  #   
  #   # print(dim(interactions.df))
  #   return(length(endpoints.targets)==0)
  # })
  
  # # Endpoints 
  # nodes_targets <- setdiff(PPI_network@interactions[,target_colname],
  #                          PPI_network@interactions[,source_colname])
  # nodes_targets <- nodes_targets[order(nodes_targets)]
  
  # non binding
  PPI_network.non_binding <- PPI_network@interactions[which(!PPI_network@interactions$arrow %in% c("---", "-+-")),]
  # setdiff(PPI_network@interactions[,target_colname], 
  #         PPI_network@interactions[,source_colname])
  nodes_targets.non_binding <- setdiff(PPI_network.non_binding[,target_colname], 
                                       PPI_network.non_binding[,source_colname])
  nodes_targets.non_binding <- nodes_targets.non_binding[order(nodes_targets.non_binding)]
  
  # remove CFTR - should be corrected with the pb TRADD
  nodes_targets.non_binding <- setdiff(nodes_targets.non_binding,
                                       "CFTR")
  # PPI_network@nodes$endpoint_topology <- PPI_network@nodes$Symbol %in% nodes_targets.non_binding
  
  
  # 
  # apoptosis_nodes <- PPI_network@nodes$Symbol[grep("CASP", PPI_network@nodes$Symbol)]
  
  # weird_endpoint.symbols <- PPI_network@nodes$Symbol[PPI_network@nodes$endpoint_topology & 
  #                                                      !PPI_network@nodes$Symbol %in% c(apoptosis_nodes, tf, "CFTR")]
  # weird_endpoint.symbols <- weird_endpoint.symbols[order(weird_endpoint.symbols)]
  
  tf <- PPI_network@nodes$Symbol[PPI_network@nodes$dorothea_tf]
  tf <- tf[order(tf)]
  
  # endpoints.binding_tf <- sapply(nodes_targets.non_binding, function(symbol){
  #   
  #   interactions.df <- PPI_network@interactions[which(PPI_network@interactions$genesymbol_source==symbol |
  #                                                       PPI_network@interactions$genesymbol_target==symbol),]
  #   binding_interactions.df <- interactions.df[which(interactions.df$arrow %in% c("---", "-+-")),]
  #   
  #   binding_nodes <- setdiff(unique(c(binding_interactions.df[,source_colname],
  #                                     binding_interactions.df[,target_colname])),
  #                            symbol)
  #   
  #   return(any(binding_nodes %in% tf))
  #   
  # })
  # 
  # 
  # endpoints.binding_source <- sapply(nodes_targets.non_binding, function(weird_endpoint){
  #   return(connected.source_or_tf(symbol = weird_endpoint,
  #                                 PPI_network = PPI_network))
  # })
  # 
  # 
  # endpoints.binding_source.2 <- sapply(nodes_targets.non_binding, function(symbol){
  #   
  #   interactions.df <- PPI_network@interactions[which(PPI_network@interactions$genesymbol_source==symbol |
  #                                                       PPI_network@interactions$genesymbol_target==symbol),]
  #   binding_interactions.df <- interactions.df[which(interactions.df$arrow %in% c("---", "-+-")),]
  #   
  #   if(dim(binding_interactions.df)[1]==0){
  #     return(FALSE)
  #   } else {
  #     binding_nodes <- unique(c(binding_interactions.df[,source_colname], 
  #                               binding_interactions.df[,target_colname]))
  #     binding_nodes <- setdiff(binding_nodes,
  #                              symbol)
  #     
  #     binding_nodes.symbols.binding_source <- sapply(binding_nodes, function(weird_endpoint){
  #       return(connected.source_or_tf(symbol = weird_endpoint,
  #                                     PPI_network = PPI_network))
  #     })
  #     return(any(binding_nodes.symbols.binding_source))
  #   }
  #   
  # })
  # 
  # nodes_targets.non_binding.final <- nodes_targets.non_binding[!(endpoints.binding_tf |
  #                                                                  endpoints.binding_source |
  #                                                                  endpoints.binding_source.2)]
  # 
  # PPI_network@nodes$endpoint_topology <- PPI_network@nodes$Symbol %in%
  #   nodes_targets.non_binding.final
  
  
  endpoints.final <- unique(c(nodes_targets.non_binding, tf))
  endpoints.final <- endpoints.final[order(endpoints.final)]
  
  PPI_network@nodes$endpoint_final <- PPI_network@nodes$Symbol %in%
    endpoints.final
  
  return(PPI_network)
  
}

# endpoints.final.df <- data.frame(Symbol = endpoints.final)
# 
# write.table(endpoints.final.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/topological_endpoint_to_phenotype_2023_04_03.tsv",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

######################################################
# TO HAVE THE ENDPOINT IN EACH PATHWAY AUTOMATICALLY #
# DON'T FORGET THAT IT IS CHECKED MANUALLY AFTERWARDS #
######################################################

# endpoints.pheno.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/topological_endpoint_to_phenotype_2023_04_03.tsv",
#                            sep = "\t",
#                            header = T,
#                            na.strings = "")
# 
# endpoints.pheno.final.df <- endpoints.pheno.df[which(endpoints.pheno.df$Endpoint_bool |
#                                                        is.na(endpoints.pheno.df$Endpoint_bool)),]
# 
# tf.final <- endpoints.pheno.final.df[which(endpoints.pheno.final.df$Endpoint_cat=="TF"),"Symbol"]
# 
# PPI_network <- CF_PPI_network.lcc.node_type
# 
# endpoints.final.pathways.df <- PPI_network@nodes[which(PPI_network@nodes$Symbol %in% endpoints.pheno.final.df$Symbol),c("Symbol", diff_pathways)]
# endpoints.final.pathways.df <- merge(endpoints.final.pathways.df,
#                                      endpoints.pheno.final.df[,c("Symbol", "Endpoint_cat")],
#                                      by = "Symbol")
# endpoints.final.pathways.df <- PPI_network@nodes[which(PPI_network@nodes$Symbol %in% tf.final),c("Symbol", diff_pathways)]
# 
# write.table(endpoints.final.pathways.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/final_endpoint_to_pathways_2023_04_04.tsv",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# endpoints.final.pathways.mat <- as.matrix(endpoints.final.pathways.df[,diff_pathways])
# rownames(endpoints.final.pathways.mat) <- endpoints.final.pathways.df$Symbol
# tf.cluster <- pheatmap(endpoints.final.pathways.mat,
#                        color = c("white", "black"),
#                        cluster_rows = F,
#                        cluster_cols = T,
#                        fontsize = 10)

######################################################
# TO HAVE THE ENDPOINT IN EACH PATHWAY AUTOMATICALLY #
# DON'T FORGET THAT IT IS CHECKED MANUALLY AFTERWARDS #
######################################################

# TF to pathways
load("/Users/matthieu/ownCloud/Thèse/Systems Biology/Transcriptomic studies/fgsea_comparison/fgsea_nes_diff_pathways_2022_09_15.RData")
diff_pathways <- rownames(fgsea_es_diff_pathways)
diff_pathways <- diff_pathways[which(diff_pathways!="AGE-RAGE signaling pathway in diabetic complications")]

endpoints.final.pathways.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/final_endpoint_to_pathways_2023_04_05_corrected.txt",
                                 sep = "\t",
                                 header = T,
                                 na.strings = "",
                                 check.names = FALSE)

endpoints.final.pathways.df[,diff_pathways] <- apply(endpoints.final.pathways.df[,diff_pathways],
                                              MARGIN = 2,
                                              function(x) as.numeric(gsub(pattern = -1,
                                                                          replacement = 1,
                                                                          x)))

# Barplot
endpoints.final.pathways.df$nb_pathways <- apply(X = endpoints.final.pathways.df[,diff_pathways], MARGIN = 1, FUN = function(x) sum(as.integer(x)))
endpoints.final.pathways.df <- endpoints.final.pathways.df[order(endpoints.final.pathways.df$nb_pathways, 
                                                                 decreasing = T),]
endpoints.final.pathways.df$Symbol <- factor(endpoints.final.pathways.df$Symbol,
                                             levels = endpoints.final.pathways.df$Symbol)

ggplot(endpoints.final.pathways.df, aes(x=Symbol, y=nb_pathways)) +
  geom_bar(stat = "identity")+
  scale_x_discrete(name = "HGNC Symbol")+
  scale_y_continuous(name = "Nb differentially expressed pathways")+
  theme(axis.text.x = element_text(face="bold", size=18, angle=90))

# Heatmap
# on peut rajouter "connected" - cf la discussion avec Véronique
endpoints.final.pathways.long <- endpoints.final.pathways.df %>%
  pivot_longer(cols = diff_pathways,
               names_to = "pathway",
               values_to = "inside")
endpoints.final.pathways.long$inside <- as.factor(endpoints.final.pathways.long$inside)
endpoints.final.pathways.heatmap <- ggplot(endpoints.final.pathways.long,aes(pathway, Symbol, fill=inside))+
  geom_tile()+
  scale_fill_manual(values = c("white", "black"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Only TF

tf.final.pathways.df <- endpoints.final.pathways.df[which(endpoints.final.df$Endpoint_cat=="TF"),]
tf.final.pathways.mat <- as.matrix(tf.final.pathways.df[,diff_pathways])
rownames(tf.final.pathways.mat) <- tf.final.pathways.df$Symbol

endpoint.df <- data.frame("cat" = tf.final.pathways.df$Endpoint_cat)
rownames(endpoint.df) <- tf.final.pathways.df$Symbol
  
tf.cluster <- pheatmap(tf.final.pathways.mat,
                       color = c("white", "black"),
                       cluster_rows = T,
                       cluster_cols = T,
                       annotation_row = endpoint.df,
                       fontsize = 10)

# # # # to test
# # source_colname="genesymbol_source"
# # target_colname="genesymbol_target"
# # PPI_network <- CF_PPI_network.lcc.curated.2.endpoint_tag
# 
# # Endpoints that are neither TF nor Capsapses
# tag_weird_endpoints <-function(PPI_network,
#                                source_colname="genesymbol_source",
#                                target_colname="genesymbol_target"){
#   
#   weird_endpoint.symbols <- PPI_network@nodes$Symbol[PPI_network@nodes$endpoint_topology & !PPI_network@nodes$dorothea_tf]
#   weird_endpoint.symbols <- setdiff(weird_endpoint.symbols,
#                                     # setdiff(weird_endpoint.symbols$Symbol,
#                                     c("CFTR"))
#   weird_endpoint.symbols <- weird_endpoint.symbols[order(weird_endpoint.symbols)]
#   
#   PPI_network@nodes$weird_endpoint <- PPI_network@nodes$Symbol %in% weird_endpoint.symbols
#   
#   weird_endpoint.interactions <- PPI_network@interactions[PPI_network@interactions$genesymbol_target %in% weird_endpoint.symbols,]
#   
#   PPI_network@interactions$weird_endpoint <- apply(X=PPI_network@interactions[,c("genesymbol_source", "genesymbol_target")],
#                                                    MARGIN = 1,
#                                                    FUN=function(x){
#                                                      return(any(x %in% PPI_network@nodes[which(PPI_network@nodes$weird_endpoint),]$Symbol))
#                                                    })
#   
#   # PPI_network@nodes$to_keep <- NULL
#   # PPI_network@interactions$to_keep <- NULL
#   
#   return(PPI_network)
#   
# }

# # Remove endpoints that are neither TF nor Capsapses
# remove_weird_endpoints <-function(PPI_network,
#                                   source_colname="genesymbol_source",
#                                   target_colname="genesymbol_target"){
#   
#   PPI_network.nodes.to_keep <- PPI_network@nodes[which(PPI_network@nodes$to_keep),]
#   
#   PPI_network.interactions.to_keep <- PPI_network@interactions[which(PPI_network@interactions$to_keep),]
#   
#   PPI_network.to_keep <- new("PPI_network",
#                              interactions=PPI_network.interactions.to_keep,
#                              nodes=PPI_network.nodes.to_keep)
#   
#   return(PPI_network)
#   
# }