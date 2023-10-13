# for get_network_nodes()
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# # to test
# PPI_network <- CF_PPI_network.lcc.node_type
# source_colname <- "genesymbol_source"
# target_colname <- "genesymbol_target"

from_PPI_to_PPathI <- function(PPI_network,
                               source_colname="genesymbol_source",
                               target_colname="genesymbol_target") {
  
  # # to test 
  # PPI_df <- interactor_to_pathway_without_borders
  # PPI_nodes <- interactor_to_pathway_without_borders_nodes
  PPI_network@nodes[which(PPI_network@nodes$CFTR_interactor),"sum"]<- 0
  PPI_network@nodes[which(PPI_network@nodes$Symbol=="CFTR"),"sum"]<- 0
  
  
  # Adding pathway occurrence to Source and Target
  PPI_network@interactions <- merge(PPI_network@interactions,
                            PPI_network@nodes[,c("Symbol", 
                                                 "sum", 
                                                 "pathway")],
                            by.x = source_colname,
                            by.y = "Symbol")
  colnames(PPI_network@interactions)[(ncol(PPI_network@interactions)-1): ncol(PPI_network@interactions)] <- 
    c("sum.source", "pathway.source")
  PPI_network@interactions <- merge(PPI_network@interactions,
                            PPI_network@nodes[,c("Symbol", 
                                                 "sum", 
                                                 "pathway")],
                            by.x = target_colname,
                            by.y = "Symbol")
  colnames(PPI_network@interactions)[(ncol(PPI_network@interactions)-1): ncol(PPI_network@interactions)] <- 
    c("sum.target", "pathway.target")
  
  # PPI_network@interactions[is.na(PPI_network@interactions$sum.source), "sum.source"] <- 0
  # PPI_network@interactions[is.na(PPI_network@interactions$sum.target), "sum.target"] <- 0
  
  col_order <- c(source_colname,
                 target_colname,
                 "effect",
                 "pathway.interaction",
                 "sum.interaction",
                 "sum.source",
                 "pathway.source",
                 "sum.target",
                 "pathway.target")
  
  PPI_interactions <- PPI_network@interactions[,col_order]
  
  intra_pathway_interactions <- apply(X = PPI_interactions[,c("sum.source", 
                                                              "sum.target")], 
                                      MARGIN = 1, 
                                      FUN = function(x) {
                                        if (any(is.na(x))){
                                          return(FALSE)
                                        } else {
                                          return(all(x==1))
                                        }
                                        })
  
  PPI_inter_pathway_interactions <- PPI_interactions[!intra_pathway_interactions,]
  
  PPathI_df <- as.data.frame(t(do.call("cbind", apply(X = PPI_inter_pathway_interactions,
                                                      MARGIN = 1,
                                                      FUN = function(x){
                                                        
                                                        print(x)
                                                        # # to test -> remove after
                                                        # x <- CFTR_to_pathway_without_borders[i_row,]
                                                        
                                                        names(x) <- colnames(PPI_inter_pathway_interactions)
                                                        output_rownames <- c("Source",
                                                                             "Target",
                                                                             "effect",
                                                                             "pathway.interaction")
                                                        
                                                        if (as.numeric(x['sum.source'])==1){
                                                          # print('hello')
                                                          final_df <-as.data.frame(c(x['pathway.source'],
                                                                                     x['genesymbol_target'],
                                                                                     x["effect"],
                                                                                     x['pathway.interaction']),
                                                                                   row.names = output_rownames)
                                                        } else if (as.numeric(x['sum.target'])==1){
                                                          # print('hi')
                                                          final_df <-as.data.frame(c(x['genesymbol_source'],
                                                                                     x['pathway.target'],
                                                                                     x["effect"],
                                                                                     x['pathway.interaction']),
                                                                                   row.names = output_rownames)
                                                        } else {
                                                          # print('coucou')
                                                          final_df <- as.data.frame(c(x["genesymbol_source"],
                                                                                      x['genesymbol_target'], 
                                                                                      x['effect'], 
                                                                                      x["pathway.interaction"]),
                                                                                    row.names = output_rownames)
                                                        }
                                                        return(final_df)
                                                      }))))
  rownames(PPathI_df) <- NULL
  
  # remove pathway --- pathway interaction
  # PPathI_df <- PPathI_df[which(PPathI_df$interaction_type.2!="intra"),]
  PPathI_df <- unique(PPathI_df)
  
  # for (i_int in which(PPathI_df$interaction_type=="CFTR_interactor_pathway")){
  #   x <- PPathI_df[i_int,]
  #   # CFTR_to_pathway_direct <- CFTR_to_pathway_without_borders[which(CFTR_to_pathway_without_borders$interaction_type=='CFTR_pathway'),]
  #   to_bind <- as.data.frame(c(x[which(x[c('Source_type','Target_type')]=="pathway")],
  #                              pathway_name,
  #                              'intra'),
  #                            col.names = c("Source",
  #                                          "Target",
  #                                          "interaction_type.2"))
  #   
  #   PPathI_df <- rbind(PPathI_df,
  #                      to_bind)
  # }
  
  PPathI_node <- get_network_nodes(network_df = PPathI_df,
                                   source_colname = "Source",
                                   target_colname = "Target")
  
  print(dim(PPathI_node))
  
  PPathI_node <- merge(PPathI_node,
                       PPI_network@nodes[,c("Symbol",
                                            "Node_type", 
                                            "CFTR_interactor", 
                                            "Type")],
                       by.x = "HGNC",
                       by.y = "Symbol",
                       all.x=TRUE)
  
  PPathI_node[which(is.na(PPathI_node$Node_type)),"Node_type"]<- "supernode"
  
  PPathI <- new("PPI_network",
                interactions=PPathI_df,
                nodes=PPathI_node)
  
  return(PPathI)
  
}

# write.table(PPathI@interactions,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_PPathI_interactions_df_2023_03_22.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(PPathI@nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_PPathI_nodes_df_2023_03_22.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)


# to test 
# PPI_network <- CF_PPI_network.scc
# test_interactions_id <- CF_PPI_network.scc@interactions[which(CF_PPI_network.scc@interactions$genesymbol_source=='LSP1' |
#                                         CF_PPI_network.scc@interactions$genesymbol_target=='LSP1'),"interaction_id"]
# test_nodes <- get_network_nodes(test_interactions,
#                                 source_colname="genesymbol_source",
#                                 target_colname="genesymbol_target")
# test_nodes <- CF_PPI_network.scc@nodes[which(CF_PPI_network.scc@nodes$Symbol %in%
#                                                test_nodes$HGNC),]
# 
# test_network <- new("PPI_network",
#                     interactions=test_interactions,
#                     nodes=test_nodes)
# PPI_network <- test_network

from_PPI_to_PPathI.only_source <- function(PPI_network,
                               source_colname="genesymbol_source",
                               target_colname="genesymbol_target") {
  
  PPI_nodes <- get_network_nodes(PPI_network@interactions,
                                 source_colname = source_colname,
                                 target_colname = target_colname)
  
  PPI_network@nodes <- merge(PPI_network@nodes,
                             PPI_nodes,
                             by.x = "Symbol",
                             by.y = "HGNC",
                             all = TRUE)
  PPI_network@nodes$Node_type <- "Protein"
  
  # Adding pathway, occurrence and degree to Source and Target
  PPI_network@interactions <- merge(PPI_network@interactions,
                                    PPI_network@nodes[,c("Symbol", 
                                                         "sum", 
                                                         "pathway",
                                                         "Count")],
                                    by.x = source_colname,
                                    by.y = "Symbol")
  colnames(PPI_network@interactions)[(ncol(PPI_network@interactions)-2): ncol(PPI_network@interactions)] <- 
    c("sum.source", "pathway.source", "degree.source")
  PPI_network@interactions <- merge(PPI_network@interactions,
                                    PPI_network@nodes[,c("Symbol", 
                                                         "sum", 
                                                         "pathway",
                                                         "Count")],
                                    by.x = target_colname,
                                    by.y = "Symbol")
  colnames(PPI_network@interactions)[(ncol(PPI_network@interactions)-2): ncol(PPI_network@interactions)] <- 
    c("sum.target", "pathway.target", "degree.target")
  
  col_order <- c(source_colname,
                 target_colname,
                 "interaction_id",
                 "effect",
                 "pathway.interaction",
                 "sum.interaction",
                 "sum.source",
                 "pathway.source",
                 "degree.source",
                 "sum.target",
                 "pathway.target",
                 "degree.target")
  
  PPI_interactions <- PPI_network@interactions[,col_order]
  
  source_nodes <- setdiff(PPI_network@interactions[,source_colname], 
                          PPI_network@interactions[,target_colname])
  # PPI_source_interactions <- PPI_interactions[which(PPI_interactions$degree.source==1),]
  
  PPathI_df <- as.data.frame(t(do.call("cbind", apply(X = PPI_interactions,
                                                      MARGIN = 1,
                                                      FUN = function(x){
                                                        
                                                        # # to test -> remove after
                                                        # x <- CFTR_to_pathway_without_borders[i_row,]
                                                        
                                                        names(x) <- colnames(PPI_interactions)
                                                        output_rownames <- c("Source",
                                                                             "Target",
                                                                             "effect",
                                                                             "pathway.interaction")
                                                        
                                                        # if (as.numeric(x['degree.source'])==1){
                                                        if (x[source_colname] %in% source_nodes & as.numeric(x['sum.source'])==1){
                                                          # print('hello')
                                                          final_df <-as.data.frame(c(x['pathway.source'],
                                                                                     x['genesymbol_target'],
                                                                                     x["effect"],
                                                                                     x['pathway.interaction']),
                                                                                   row.names = output_rownames)
                                                        # } else if (as.numeric(x['sum.target'])==1){
                                                        #   # print('hi')
                                                        #   final_df <-as.data.frame(c(x['genesymbol_source'],
                                                        #                              x['pathway.target'],
                                                        #                              x["effect"],
                                                        #                              x['pathway.interaction']),
                                                        #                            row.names = output_rownames)
                                                        } else {
                                                          # print('coucou')
                                                          final_df <- as.data.frame(c(x["genesymbol_source"],
                                                                                      x['genesymbol_target'], 
                                                                                      x['effect'], 
                                                                                      x["pathway.interaction"]),
                                                                                    row.names = output_rownames)
                                                        }
                                                        return(final_df)
                                                      }))))
  rownames(PPathI_df) <- NULL
  
  # remove pathway --- pathway interaction
  # PPathI_df <- PPathI_df[which(PPathI_df$interaction_type.2!="intra"),]
  PPathI_df <- unique(PPathI_df)
  
  # for (i_int in which(PPathI_df$interaction_type=="CFTR_interactor_pathway")){
  #   x <- PPathI_df[i_int,]
  #   # CFTR_to_pathway_direct <- CFTR_to_pathway_without_borders[which(CFTR_to_pathway_without_borders$interaction_type=='CFTR_pathway'),]
  #   to_bind <- as.data.frame(c(x[which(x[c('Source_type','Target_type')]=="pathway")],
  #                              pathway_name,
  #                              'intra'),
  #                            col.names = c("Source",
  #                                          "Target",
  #                                          "interaction_type.2"))
  #   
  #   PPathI_df <- rbind(PPathI_df,
  #                      to_bind)
  # }
  
  PPathI_node <- get_network_nodes(network_df = PPathI_df,
                                   source_colname = "Source",
                                   target_colname = "Target")
  
  print(dim(PPathI_node))
  
  PPathI_node <- merge(PPathI_node,
                       PPI_network@nodes[,c("Symbol","Node_type", "sum", "pathway", "endpoint_to_keep")],
                       by.x = "HGNC",
                       by.y = "Symbol",
                       all.x=TRUE)
  
  PPathI_node[which(is.na(PPathI_node$Node_type)),"Node_type"]<- "supernode"
  PPathI_node[which(PPathI_node$Node_type=="supernode"), "pathway"] <- PPathI_node[which(PPathI_node$Node_type=="supernode"), "HGNC"]
  PPathI_node[which(PPathI_node$Node_type=="supernode"), "sum"] <- 1
  PPathI_node[which(PPathI_node$Node_type=="supernode"), "endpoint_to_keep"] <- FALSE
  PPathI <- new("PPI_network",
                interactions=PPathI_df,
                nodes=PPathI_node)
  
  return(PPathI)
  
}
