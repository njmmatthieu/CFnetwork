# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

library(dplyr)
library(ggplot2)
library(igraph)
library(pheatmap)

load("/Users/matthieu/ownCloud/Thèse/Systems Biology/Transcriptomic studies/fgsea_comparison/fgsea_nes_diff_pathways_2022_09_15.RData")
diff_pathways <- rownames(fgsea_es_diff_pathways)
diff_pathways <- diff_pathways[which(diff_pathways!="AGE-RAGE signaling pathway in diabetic complications")]

# # to test
# protein <- "HSP90AA1"

CF_PPI_network.lcc.node_type.interactions <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_interactions_df_2023_04_18.txt",
            sep = "\t",
            header = T,
            check.names = F)

CF_PPI_network.lcc.node_type.nodes <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_nodes_df_2022_04_18.txt",
            sep = "\t",
            header = T,
            check.names = F)

CF_PPI_network.lcc.node_type <- new("PPI_network",
                      interactions=CF_PPI_network.lcc.node_type.interactions,
                      nodes=CF_PPI_network.lcc.node_type.nodes)

# Remove CFTR from the network
CF_PPI_network.lcc.node_type@interactions <- CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$genesymbol_source!="CFTR" &
                                                                                               CF_PPI_network.lcc.node_type@interactions$genesymbol_target!="CFTR"),]
CF_PPI_network.lcc.node_type@nodes <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes$Symbol!="CFTR"),]

##################
# MODIFIER GENES #
##################

genes_in_litterature <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/CF_genes_in_literature/CF_genes_in_literature_2v_HGNC_2021_10_28.csv",
                                   sep = "\t",
                                   header = T)

#############
# FUNCTIONS #
#############

# # to test
# network <- CF_PPI_network.llc.curated.2.node_type
# protein <- "TRADD"

# Get the interactions of the protein in the network
interactions_of_protein <- function(protein, 
                                   network) {
  
  interactions.source <- network@interactions[which(network@interactions$genesymbol_source==protein),]
  
  interactions.target <- network@interactions[which(network@interactions$genesymbol_target==protein),]
  
  interactions.df <- rbind(interactions.source,
                           interactions.target)
  
  interactions.unique.df <- unique(interactions.df)
  
  return(interactions.unique.df)
  
}



# Get the pathways connected to the protein
pathways_connected_function <- function(protein, network) {
  
  # Targets of the protein
  protein.targets_interactions <- network@interactions[which(network@interactions$genesymbol_source==protein),c("effect",
                                                                                                   "arrow",
                                                                                                   "genesymbol_source",
                                                                                                   "genesymbol_target")]
  # Remove CFTR interactors targets
  CFTR_interactors <- network@nodes[which(network@nodes$CFTR_interactor),"Symbol"]
  protein.targets_interactions.non_CFTR_interactors <-   protein.targets_interactions[which(!protein.targets_interactions$genesymbol_target %in% CFTR_interactors),]
  targets.non_CFTR_interactors <- protein.targets_interactions.non_CFTR_interactors$genesymbol_target
  # print(protein.targets_interactions)
  
  # Retrieving pathways of the targets
  protein.targets <- network@nodes[which(network@nodes$Symbol %in% targets.non_CFTR_interactors),c("Symbol",
                                                                                                   diff_pathways,
                                                                                                   "pathway")]
  
  protein.targets.df <- merge(protein.targets_interactions.non_CFTR_interactors,
                              protein.targets,
                              by.x = "genesymbol_target",
                              by.y = "Symbol")
  
  protein.targets.df <- protein.targets.df[,c("genesymbol_source",
                                              "genesymbol_target",
                                              "effect",
                                              "arrow",
                                              diff_pathways)]
                                              # "pathway")]
  
  protein.targets.df[,diff_pathways] <- sapply(protein.targets.df[,diff_pathways], function(x){
    x[is.na(x)] <- 0
    return(x)
  })
  
  # Wide to long format
  protein.targets.df.long <- pivot_longer(data = protein.targets.df, cols = diff_pathways, names_to = "diff_pathway", values_to = "occurrence")
  protein.targets.df.long <- protein.targets.df.long[which(protein.targets.df.long$occurrence==1),c("genesymbol_source", 
                                                                                                    "diff_pathway", 
                                                                                                    "effect", 
                                                                                                    "arrow")]
  # peut-être changer le nom des colonnes
  colnames(protein.targets.df.long) <- c("genesymbol_source", 
                                         "genesymbol_target", 
                                         "effect", 
                                         "arrow")
  # print(protein.targets.df)
  # Sources of the protein
  protein.sources_interactions <- network@interactions[which(network@interactions$genesymbol_target==protein),c("effect",
                                                                                                   "arrow",
                                                                                                   "genesymbol_source",
                                                                                                   "genesymbol_target")]
  
  # print(protein.sources_interactions)
  protein.sources_interactions.non_CFTR_interactors <-   protein.sources_interactions[which(!protein.sources_interactions$genesymbol_source %in% CFTR_interactors),]
  sources.non_CFTR_interactors <- protein.sources_interactions.non_CFTR_interactors$genesymbol_source
  # print(protein.sources_interactions.non_CFTR_interactors)
  
  protein.sources <- network@nodes[which(network@nodes$Symbol %in% sources.non_CFTR_interactors),c("Symbol",
                                                                                                   diff_pathways,
                                                                                                   "pathway")]
  
  protein.sources.df <- merge(protein.sources_interactions.non_CFTR_interactors,
                              protein.sources,
                              by.x = "genesymbol_source",
                              by.y = "Symbol")
  
  protein.sources.df <- protein.sources.df[,c("genesymbol_source",
                                              "genesymbol_target",
                                              "effect",
                                              "arrow",
                                              diff_pathways)]
                                              # "pathway")]
  print(protein.sources.df)
  
  # print(protein.sources.df)
  protein.sources.df[,diff_pathways] <- sapply(protein.sources.df[,diff_pathways], function(x){
    x[is.na(x)] <- 0
    return(x)
  })
  
  protein.sources.df.long <- pivot_longer(data = protein.sources.df, 
                                          cols = diff_pathways, 
                                          names_to = "diff_pathway", 
                                          values_to = "occurrence")
  protein.sources.df.long <- protein.sources.df.long[which(protein.sources.df.long$occurrence==1),c("diff_pathway",
                                                                                                    "genesymbol_target",
                                                                                                    "effect", 
                                                                                                    "arrow")]
  
  # peut-être changer le nom des colonnes
  colnames(protein.sources.df.long) <- c("genesymbol_source", 
                                         "genesymbol_target", 
                                         "effect", 
                                         "arrow")
  
  
  protein.interactors <- rbind(protein.targets.df.long,
                               protein.sources.df.long)
  
  return(unique(protein.interactors))
  
}

test <- pathways_connected_function("SLC9A3R1",
                                    CF_PPI_network.lcc.node_type)
# 
# test.2 <- pathways_connected_function("STUB1",
#                                    CF_PPI_network.llc.curated.2.node_type)

#######################################
# CFTR interactors - pathways #
#######################################


# 1 - CFTR interactor

CFTR_interactor.in_CF_PPI_network <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes[,"CFTR_interactor"]),]

all_CFTR_interactors.pathways_connected <- lapply(CFTR_interactor.in_CF_PPI_network$Symbol, function(CFTR.interactor){
  return(pathways_connected_function(CFTR.interactor, 
                                     CF_PPI_network.lcc.node_type))
})

all_CFTR_interactors.pathways_connected.df <- do.call("rbind", all_CFTR_interactors.pathways_connected)

# 2 - CFTR interactors - CFTR interactors interactions

CFTR_interactors.interactions <- CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$genesymbol_source %in% CFTR_interactor.in_CF_PPI_network$Symbol &
                                                                                   CF_PPI_network.lcc.node_type@interactions$genesymbol_target %in% CFTR_interactor.in_CF_PPI_network$Symbol),]

all_CFTR_interactors.pathways_connected.interactions <- rbind(all_CFTR_interactors.pathways_connected.df,
                                                       CFTR_interactors.interactions[,c("genesymbol_source",
                                                                                        "genesymbol_target",
                                                                                        "effect",
                                                                                        "arrow")])

# Node type

all_CFTR_interactors.pathways_connected.nodes <- get_network_nodes(network_df = data.frame(all_CFTR_interactors.pathways_connected.interactions),
                                 source_colname = "genesymbol_source",
                                 target_colname = "genesymbol_target")

all_CFTR_interactors.pathways_connected.nodes <- merge(all_CFTR_interactors.pathways_connected.nodes,
                                                       CF_PPI_network.lcc.node_type@nodes[,c("Symbol",
                                          "Node_type", 
                                          "CFTR_interactor", 
                                          "status")],
                     by.x = "HGNC",
                     by.y = "Symbol",
                     all.x=TRUE)

all_CFTR_interactors.pathways_connected.nodes[which(is.na(all_CFTR_interactors.pathways_connected.nodes$Node_type)),"Node_type"]<- "supernode"


# write.table(all_CFTR_interactors.pathways_connected.interactions,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_PPathI_interactions_df_2023_04_07.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(all_CFTR_interactors.pathways_connected.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_PPathI_nodes_df_2023_04_07.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

CFTR_interactors_connected_to_kegg_diff_pathways <- data.frame(sapply(CFTR_interactor.in_CF_PPI_network$Symbol, function(interactor){
  
  sapply(diff_pathways, function(pathway_name){
    
    print(pathway_name)
    # # to test
    pathway_genes <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes[,pathway_name]==1),"Symbol"]
    print(length(pathway_genes))
    
    if (interactor %in% pathway_genes){
      return("inside - pathway")
    } else {
      
      interactor.interactions <- interactions_of_protein(protein = interactor,
                                                         network = CF_PPI_network.lcc.node_type)
      
      if (any(interactor.interactions$genesymbol_target %in% pathway_genes) |
          any(interactor.interactions$genesymbol_source %in% pathway_genes)){
        return("connected - pathway")
      } else {
        return("not connected")
      }
      
    }
  })
}) 
)

# intesting_interactors <- colnames(CFTR_interactors_connected_to_kegg_diff_pathways)
CFTR_interactors_connected_to_kegg_diff_pathways$pathway_name <- rownames(CFTR_interactors_connected_to_kegg_diff_pathways)
CFTR_interactors_connected_to_kegg_diff_pathways.long <- pivot_longer(CFTR_interactors_connected_to_kegg_diff_pathways,
                                                                      cols = CFTR_interactor.in_CF_PPI_network$Symbol,
                                                                      names_to = "Interactor",
                                                                      values_to = "state")
CFTR_interactors_connected_to_kegg_diff_pathways.long$state <- factor(CFTR_interactors_connected_to_kegg_diff_pathways.long$state,
                                                                      levels = c("inside - pathway",
                                                                                    "connected - pathway",
                                                                                     "not connected"))
CFTR_interactors_connected_to_kegg_diff_pathways.long$state_value <- as.integer(CFTR_interactors_connected_to_kegg_diff_pathways.long$state)-1

# Heatmap

CFTR_interactors_connected_to_kegg_diff_pathways.heatmap <- 
  ggplot(CFTR_interactors_connected_to_kegg_diff_pathways.long, 
         aes(y = pathway_name, 
             x = Interactor,
             fill = state))+
  geom_tile(color = "white",
            lwd = 1.5,
            linetype = 1)+
  scale_fill_manual(values = c("#FFB570","#FFD966","#D9D9D9"))+
  # values = c("inside","connected","not connected"))+
  theme(panel.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10, vjust =1, hjust = 1, angle = 45),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text(size=24))

# Pheatmap

CFTR_interactors_connected_to_kegg_diff_pathways.int <- data.frame(sapply(CFTR_interactors_connected_to_kegg_diff_pathways, function(x) {
  return(as.integer(factor(x, levels= rev(c("inside - pathway",
                                            "connected - pathway",
                                            "not connected"))))-1)
}))
rownames(CFTR_interactors_connected_to_kegg_diff_pathways.int) <-
  rownames(CFTR_interactors_connected_to_kegg_diff_pathways)
pheatmap(CFTR_interactors_connected_to_kegg_diff_pathways.int,
         cluster_cols = T,
         cluster_rows = T)

#############################################
# CFTR interactors - intermediate - pathway #
#############################################

CFTR_interactor.in_CF_PPI_network <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes[,"CFTR_interactor"]),]

CFTR_interactors.intermediate.interactions <- CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$genesymbol_source %in% CFTR_interactor.in_CF_PPI_network$Symbol |
                                                                                                CF_PPI_network.lcc.node_type@interactions$genesymbol_target %in% CFTR_interactor.in_CF_PPI_network$Symbol),]

# Find CFTR interactors intermediates
CFTR_interactors.intermediates <- unique(c(CFTR_interactors.intermediate.interactions$genesymbol_source,
                                         CFTR_interactors.intermediate.interactions$genesymbol_target))

CFTR_interactors.intermediates <- setdiff(CFTR_interactors.intermediates,
                                          CFTR_interactor.in_CF_PPI_network$Symbol)

# # Directed

# all_CFTR_interactors.intermediates.pathways_connected <- lapply(CFTR_interactors.intermediates, function(intermediate){
#   return(pathways_connected_function(intermediate, 
#                                      CF_PPI_network.llc.curated.2.node_type))
# })
# all_CFTR_interactors.intermediates.pathways_connected.df <- do.call("rbind", all_CFTR_interactors.intermediates.pathways_connected)
# 

# Undirected
all_CFTR_interactors.intermediates.pathways_connected.wide <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes$Symbol %in% CFTR_interactors.intermediates),c("Symbol",
                                                                                                                                              diff_pathways)]
all_CFTR_interactors.intermediates.pathways_connected.df <- pivot_longer(data = all_CFTR_interactors.intermediates.pathways_connected.wide,
                                                                           cols = diff_pathways,
                                                                           names_to = "genesymbol_target",
                                                                           values_to = "inside")

# Heatmap
rownames(all_CFTR_interactors.intermediates.pathways_connected.wide) <- all_CFTR_interactors.intermediates.pathways_connected.wide$Symbol
all_CFTR_interactors.intermediates.pathways_connected.wide$Symbol <- NULL

pheatmap(all_CFTR_interactors.intermediates.pathways_connected.wide,
         cluster_cols = T,
         cluster_rows = T)


all_CFTR_interactors.intermediates.pathways_connected.df <- all_CFTR_interactors.intermediates.pathways_connected.df[which(all_CFTR_interactors.intermediates.pathways_connected.df$inside!=0),]
all_CFTR_interactors.intermediates.pathways_connected.df$inside <- NULL
colnames(all_CFTR_interactors.intermediates.pathways_connected.df)[1] <- "genesymbol_source"
all_CFTR_interactors.intermediates.pathways_connected.df$effect <- "undirected"
all_CFTR_interactors.intermediates.pathways_connected.df$arrow <- "NA"

all_CFTR_interactors.intermediates.pathways_connected.interactions <- rbind(all_CFTR_interactors.intermediates.pathways_connected.df,
                                                                            CFTR_interactors.intermediate.interactions[,c("genesymbol_source",
                                                                                               "genesymbol_target",
                                                                                               "effect",
                                                                                               "arrow")])
all_CFTR_interactors.intermediates.pathways_connected.nodes <- get_network_nodes(network_df = data.frame(all_CFTR_interactors.intermediates.pathways_connected.interactions),
                                                                   source_colname = "genesymbol_source",
                                                                   target_colname = "genesymbol_target")

all_CFTR_interactors.intermediates.pathways_connected.nodes <- merge(all_CFTR_interactors.intermediates.pathways_connected.nodes,
                                                                     CF_PPI_network.lcc.node_type@nodes[,c("Symbol",
                                                                                                       "Node_type", 
                                                                                                       "CFTR_interactor", 
                                                                                                       "status")],
                                                       by.x = "HGNC",
                                                       by.y = "Symbol",
                                                       all.x=TRUE)

all_CFTR_interactors.intermediates.pathways_connected.nodes[which(is.na(all_CFTR_interactors.intermediates.pathways_connected.nodes$Node_type)),"Node_type"]<- "supernode"

# write.table(all_CFTR_interactors.intermediates.pathways_connected.interactions,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_and_intermediates_PPathI_interactions_undirected_2023_04_07.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(all_CFTR_interactors.intermediates.pathways_connected.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_and_intermediates_PPathI_nodes_undirected_2023_04_07.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)




# # Undirected all nodes
# all_nodes.pathways_connected.df <- pivot_longer(data = CF_PPI_network.llc.curated.2.node_type@nodes[,c("Symbol",
#                                                                                                                                 diff_pathways)],
#                                                                          cols = diff_pathways,
#                                                                          names_to = "genesymbol_target",
#                                                                          values_to = "inside")
# all_nodes.pathways_connected.df <- all_nodes.pathways_connected.df[which(all_nodes.pathways_connected.df$inside!=0),]
# all_nodes.pathways_connected.df$inside <- NULL
# colnames(all_nodes.pathways_connected.df)[1] <- "genesymbol_source"
# all_nodes.pathways_connected.df$effect <- "undirected"
# all_nodes.pathways_connected.df$arrow <- "NA"
# 
# all_nodes.pathways_connected.nodes <- get_network_nodes(network_df = data.frame(all_nodes.pathways_connected.df),
#                                                                                  source_colname = "genesymbol_source",
#                                                                                  target_colname = "genesymbol_target")
# 
# all_nodes.pathways_connected.nodes <- merge(all_nodes.pathways_connected.nodes,
#                                                                      CF_PPI_network.llc.curated.2.node_type@nodes[,c("Symbol",
#                                                                                                                      "Node_type", 
#                                                                                                                      "CFTR_interactor", 
#                                                                                                                      "Type")],
#                                                                      by.x = "HGNC",
#                                                                      by.y = "Symbol",
#                                                                      all.x=TRUE)
# 
# all_nodes.pathways_connected.nodes[which(is.na(all_nodes.pathways_connected.nodes$Node_type)),"Node_type"]<- "supernode"
# 
# write.table(all_nodes.pathways_connected.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/all_nodes_PPathI_interactions_undirected_2023_03_29.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(all_nodes.pathways_connected.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/all_nodes_PPathI_nodes_undirected_2023_03_29.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

# Undirected
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

endpoits.pathways_connected.df <- pivot_longer(data = endpoints.final.pathways.df,
                                               cols = diff_pathways,
                                               names_to = "genesymbol_target",
                                               values_to = "inside")

endpoits.pathways_connected.df <- endpoits.pathways_connected.df[which(endpoits.pathways_connected.df$inside!=0),]

endpoits.pathways_connected.df$inside <- NULL
colnames(endpoits.pathways_connected.df)[1] <- "genesymbol_source"
endpoits.pathways_connected.df$effect <- "undirected"
endpoits.pathways_connected.df$arrow <- "NA"

endpoits.pathways_connected.nodes <- get_network_nodes(network_df = data.frame(endpoits.pathways_connected.df),
                                                                                 source_colname = "genesymbol_source",
                                                                                 target_colname = "genesymbol_target")

endpoits.pathways_connected.nodes <- merge(endpoits.pathways_connected.nodes,
                                           CF_PPI_network.lcc.node_type@nodes[,c("Symbol",
                                                                                 "Node_type", 
                                                                                 "CFTR_interactor", 
                                                                                 "status")],
                                           by.x = "HGNC",
                                           by.y = "Symbol",
                                           all.x=TRUE)

endpoits.pathways_connected.nodes[which(is.na(endpoits.pathways_connected.nodes$Node_type)),"Node_type"]<- "supernode"

# write.table(endpoits.pathways_connected.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/endpoints_PPathI_interactions_undirected_2023_04_10.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(endpoits.pathways_connected.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/endpoints_PPathI_nodes_undirected_2023_04_10.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)


# Most connected nodes

CF_PPI_network.lcc.node_type.nodes.degree <- get_network_nodes(CF_PPI_network.lcc.node_type@interactions,
                                                               source_colname = "genesymbol_source",
                                                               target_colname = "genesymbol_target")
CF_PPI_network.lcc.node_type.nodes.degree <- CF_PPI_network.lcc.node_type.nodes.degree[order(CF_PPI_network.lcc.node_type.nodes.degree$Count,
                                                                                             decreasing = T),]

# Basic histogram
ggplot(CF_PPI_network.lcc.node_type.nodes.degree, aes(x=Count)) +
  geom_histogram(binwidth=1)+
  scale_x_continuous(name = "Degree")+
  scale_y_continuous(name = "Nb of proteins")

# 30 first1
CF_PPI_network.lcc.node_type.nodes.degree.9more <- 
  CF_PPI_network.lcc.node_type.nodes.degree %>%
  filter(Count>=9)
CF_PPI_network.lcc.node_type.nodes.degree.9more$HGNC <-
  factor(CF_PPI_network.lcc.node_type.nodes.degree.9more$HGNC,
         levels = CF_PPI_network.lcc.node_type.nodes.degree.9more$HGNC)

# Barplot
ggplot(CF_PPI_network.lcc.node_type.nodes.degree.9more, aes(x=HGNC, y=Count)) +
  geom_bar(stat = "identity")+
  scale_x_discrete(name = "HGNC Symbol")+
  scale_y_continuous(name = "Degree")+
  theme(axis.text.x = element_text(face="bold", size=18, angle=90))

# Heatmap

CF_PPI_network.nodes.degree.9more.pathways_connected.wide <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes$Symbol %in% CF_PPI_network.lcc.node_type.nodes.degree.9more$HGNC),c("Symbol",diff_pathways)]
CF_PPI_network.nodes.degree.9more.pathways_connected.wide$Symbol <-
  factor(CF_PPI_network.nodes.degree.9more.pathways_connected.wide$Symbol,
         levels = CF_PPI_network.lcc.node_type.nodes.degree.9more$HGNC)
CF_PPI_network.nodes.degree.9more.pathways_connected.wide <- t(CF_PPI_network.nodes.degree.9more.pathways_connected.wide[order(CF_PPI_network.nodes.degree.9more.pathways_connected.wide$Symbol),])

# After correction

hubs.final.pathways.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/hubs_to_pathways_2023_04_13_corrected.txt",
                                          sep = "\t",
                                          header = T,
                                          na.strings = "",
                                          check.names = FALSE)
hubs.final.pathways.df[,diff_pathways] <- apply(hubs.final.pathways.df[,diff_pathways],
                                              MARGIN = 2,
                                              function(x) as.numeric(gsub(pattern = -1,
                                                                          replacement = 1,
                                                                          x)))

# Heatmap

hubs.final.pathways.mat <- as.matrix(hubs.final.pathways.df[,diff_pathways])
rownames(hubs.final.pathways.mat) <- hubs.final.pathways.df$Symbol

# endpoint.df <- data.frame("cat" = tf.final.pathways.df$Endpoint_cat)
# rownames(endpoint.df) <- tf.final.pathways.df$Symbol

hubs.cluster <- pheatmap(t(hubs.final.pathways.mat),
                       color = c("white", "black"),
                       cluster_rows = F,
                       cluster_cols = F,
                       show_colnames = F,
                       fontsize = 10)


hubs.final.pathways.long <- pivot_longer(data = hubs.final.pathways.df,
                                               cols = diff_pathways,
                                               names_to = "genesymbol_target",
                                               values_to = "inside")

hubs.final.pathways.long <- hubs.final.pathways.long[which(hubs.final.pathways.long$inside!=0),]

hubs.final.pathways.long$inside <- NULL
colnames(hubs.final.pathways.long)[1] <- "genesymbol_source"
hubs.final.pathways.long$effect <- "undirected"
hubs.final.pathways.long$arrow <- "NA"

hubs.interacctions.bool <- apply(CF_PPI_network.lcc.node_type@interactions[,c("genesymbol_source", "genesymbol_target")], MARGIN=1, FUN = function(x){
  return(all(x %in% hubs.final.pathways.df$Symbol))
})
hubs.interactions.df <- CF_PPI_network.lcc.node_type@interactions[hubs.interacctions.bool,]
# write.table(hubs.interactions.df,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/hubs_interactions_2023_04_14.tsv",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

hubs.interactions.PPathI <- rbind(hubs.interactions.df[,c("genesymbol_source",
                                                          "genesymbol_target",
                                                          "effect",
                                                          "arrow")], 
                                  hubs.final.pathways.long)

hubs.interactions.PPathI.nodes <- get_network_nodes(network_df = hubs.interactions.PPathI,
                                                       source_colname = "genesymbol_source",
                                                       target_colname = "genesymbol_target")

hubs.interactions.PPathI.nodes <- merge(hubs.interactions.PPathI.nodes,
                                           CF_PPI_network.lcc.node_type@nodes[,c("Symbol",
                                                                                 "Node_type", 
                                                                                 "CFTR_interactor", 
                                                                                 "status")],
                                           by.x = "HGNC",
                                           by.y = "Symbol",
                                           all.x=TRUE)

hubs.interactions.PPathI.nodes[which(is.na(hubs.interactions.PPathI.nodes$Node_type)),"Node_type"]<- "supernode"


# write.table(hubs.interactions.PPathI,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/hubs_PPathI_interactions_undirected_2023_04_14.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(hubs.interactions.PPathI.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/hubs_PPathI_nodes_undirected_2023_04_14.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

hubs.interactions.PPathI.for_igraph <- hubs.interactions.PPathI[,c("genesymbol_source",
                                                                   "genesymbol_target")]
colnames(hubs.interactions.PPathI.for_igraph) <- c("from",
                                                   "to")

hubs.interactions.PPathI.igraph <- graph_from_data_frame(hubs.interactions.PPathI.for_igraph, directed=FALSE)

hubs.interactions.PPathI.bc <- as.data.frame(betweenness(hubs.interactions.PPathI.igraph))

# All network

library(tidyverse)
library(dplyr)

hubs.final.pathways.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/hubs_to_pathways_2023_04_13_corrected.txt",
                                     sep = "\t",
                                     header = T,
                                     na.strings = "",
                                     check.names = FALSE)

endpoints.final.pathways.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/final_endpoint_to_pathways_2023_04_05_corrected.txt",
                                          sep = "\t",
                                          header = T,
                                          na.strings = "",
                                          check.names = FALSE)



## Interesting nodes
hubs <- as.character(hubs.final.pathways.df$Symbol)
endpoints.final.pathways.df <- endpoints.final.pathways.df[order(endpoints.final.pathways.df$Endpoint_cat),]
endpoints.final.pathways.df$Symbol <- as.character(endpoints.final.pathways.df$Symbol)
endpoints <- as.character(endpoints.final.pathways.df$Symbol)
tf <- CF_PPI_network.lcc.node_type@nodes$Symbol[CF_PPI_network.lcc.node_type@nodes$dorothea_tf]
tf <- tf[order(tf)]
# CFTR_interactors <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes$CFTR_interactor),"Symbol"]
CFTR_interactors <- c("HSP90AA1",
                      "TRADD",
                      "PRKACA",
                      "SYK",
                      "CSNK2A1",
                      "SRC",
                      "PLCB1",
                      "EZR")

endpoint.cat.annot <- data.frame("cat" = endpoints.final.pathways.df$Endpoint_cat)
rownames(endpoint.cat.annot) <- endpoints.final.pathways.df$Symbol

## Pb of binding interactions
## Non binding
CF_PPI_network.lcc.node_type.interactions.non_binding <- CF_PPI_network.lcc.node_type@interactions[which(!(CF_PPI_network.lcc.node_type@interactions$effect %in% c("binding/association"))),]
CF_PPI_network.lcc.node_type.interactions.non_binding <- CF_PPI_network.lcc.node_type.interactions.non_binding[,c("genesymbol_source",
                                                                                                                  "genesymbol_target")]
colnames(CF_PPI_network.lcc.node_type.interactions.non_binding) <- c("from", "to")

# ## Dissociation
# CF_PPI_network.lcc.node_type.interactions.dissociation <- CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$effect %in% c("dissociation")),]

CF_PPI_network.lcc.node_type.interactions.binding <- CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$effect %in% c("binding/association")),]
# Both directions for binding interactions
## one direction
CF_PPI_network.lcc.node_type.interactions.binding.one_direction <- CF_PPI_network.lcc.node_type.interactions.binding[,c("genesymbol_source","genesymbol_target")]
colnames(CF_PPI_network.lcc.node_type.interactions.binding.one_direction) <- c("from", "to")
## other direction
CF_PPI_network.lcc.node_type.interactions.binding.other_direction <- CF_PPI_network.lcc.node_type.interactions.binding[,c("genesymbol_target", "genesymbol_source")]
colnames(CF_PPI_network.lcc.node_type.interactions.binding.other_direction) <- c("from", "to")
## both directions
CF_PPI_network.lcc.node_type.interactions.binding.both_directions <- rbind(CF_PPI_network.lcc.node_type.interactions.binding.one_direction,
                                                                           CF_PPI_network.lcc.node_type.interactions.binding.other_direction)


# # check for TF binding non-TF
# test <- CF_PPI_network.lcc.node_type.interactions.binding.both_directions[apply(X = CF_PPI_network.lcc.node_type.interactions.binding.both_directions,MARGIN = 1, FUN = function(x) sum(x %in% tf)==1),]


CF_PPI_network.lcc.for_igraph <- rbind(CF_PPI_network.lcc.node_type.interactions.non_binding,
                                       CF_PPI_network.lcc.node_type.interactions.binding.both_directions)

# Pb of CFTR interactors
CFTR_interactors.Symbol <- CF_PPI_network.lcc.node_type@nodes$Symbol[which(CF_PPI_network.lcc.node_type@nodes$CFTR_interactor)]
CFTR_interactors.non_source <- CFTR_interactors.Symbol[!CFTR_interactors.Symbol %in% CF_PPI_network.lcc.for_igraph$from]
CFTR_interactors.non_source.for_igraph <- CF_PPI_network.lcc.for_igraph[which(CF_PPI_network.lcc.for_igraph$to %in% CFTR_interactors.non_source),]
CFTR_interactors.non_source.for_igraph.other_direction <- CFTR_interactors.non_source.for_igraph[,c("to",
                                                                                                    "from")]
colnames(CFTR_interactors.non_source.for_igraph.other_direction)<- c("from",
                                                                     "to")
CF_PPI_network.lcc.for_igraph.2 <- rbind(CF_PPI_network.lcc.for_igraph,
                                         CFTR_interactors.non_source.for_igraph.other_direction)

# Igraph
CF_PPI_network.lcc.igraph <- graph_from_data_frame(CF_PPI_network.lcc.for_igraph.2, directed=TRUE)

## Distances
CF_PPI_network.lcc.dists <- distances(CF_PPI_network.lcc.igraph, mode = "out")

### CFTR interactors to Endpoints
CF_PPI_network.lcc.interactors_to_endpoints <- CF_PPI_network.lcc.dists[CFTR_interactors, endpoints]

#### remove CSNK2A1
CF_PPI_network.lcc.interactors_to_endpoints.2 <- CF_PPI_network.lcc.interactors_to_endpoints[which(row.names(CF_PPI_network.lcc.interactors_to_endpoints)!="CSNK2A1"), ]

##### Heatmap
CF_PPI_network.lcc.interactors_to_endpoints.heatmap <- pheatmap(CF_PPI_network.lcc.interactors_to_endpoints.2,
                         # color = c("white", "black"),
                         cluster_rows = T,
                         cluster_cols = F,
                         # show_colnames = T,
                         annotation_col = endpoint.cat.annot,
                         fontsize = 10)
##### Barplot
CF_PPI_network.lcc.interactors_to_endpoints.df <- as.data.frame(CF_PPI_network.lcc.interactors_to_endpoints)
CF_PPI_network.lcc.interactors_to_endpoints.df$CFTR_interactor <- rownames(CF_PPI_network.lcc.interactors_to_endpoints.df)
CF_PPI_network.lcc.interactors_to_endpoints.df.longer <- pivot_longer(CF_PPI_network.lcc.interactors_to_endpoints.df,
                                                                     cols = endpoints,
                                                                     names_to = "endpoint",
                                                                     values_to = "distance")
CF_PPI_network.lcc.interactors_to_endpoints.df.longer <- merge(CF_PPI_network.lcc.interactors_to_endpoints.df.longer,
                                                               endpoints.final.pathways.df[,c("Symbol", "Endpoint_cat")],
                                                               by.x = "endpoint",
                                                               by.y = "Symbol",
                                                               all.x = TRUE)

CF_PPI_network.lcc.interactors_to_endpoints.barplot <- ggplot(CF_PPI_network.lcc.interactors_to_endpoints.df.longer, aes(x=distance, fill = Endpoint_cat))+
  geom_bar(stat="bin")+
  facet_wrap(~ CFTR_interactor)

# only Regulation of actin cytoskeleton

CF_PPI_network.lcc.interactors_to_REgActin_endpoints.df.longer <- CF_PPI_network.lcc.interactors_to_endpoints.df.longer[which(CF_PPI_network.lcc.interactors_to_endpoints.df.longer$Endpoint_cat=="Regulation of actin cytoskeleton"),]
CF_PPI_network.lcc.interactors_to_REgActin_endpoints.barplot <- ggplot(CF_PPI_network.lcc.interactors_to_REgActin_endpoints.df.longer, aes(x=distance))+
  geom_bar(stat="bin")+
  facet_wrap(~ CFTR_interactor)

# Mean by interactor

CF_PPI_network.lcc.interactors_to_endpoints.df.longer %>% 
  group_by(CFTR_interactor) %>%
  summarise(across(distance, mean, na.rm=TRUE))


CF_PPI_network.lcc.interactors_to_endpoints.2 <- CF_PPI_network.lcc.dists[endpoints, CFTR_interactors]

which(names(V(CF_PPI_network.lcc.igraph))=="CSNK2A1")
# 144
which(names(V(CF_PPI_network.lcc.igraph))=="TRADD")
# 87
which(names(V(CF_PPI_network.lcc.igraph))=="SYK")
# 329
which(names(V(CF_PPI_network.lcc.igraph))=="ACTN4")
# 312
which(names(V(CF_PPI_network.lcc.igraph))=="ARPC5")
# 230
which(names(V(CF_PPI_network.lcc.igraph))=="EZR")

print(all_shortest_paths(CF_PPI_network.lcc.igraph,
                         from=V(CF_PPI_network.lcc.igraph)[230],
                         to=V(CF_PPI_network.lcc.igraph)[312],
                         mode="out")$res)


### CFTR interactors to Hubs
CF_PPI_network.lcc.interactors_to_hubs <- CF_PPI_network.lcc.dists[CFTR_interactors, hubs]
#### remove CSNK2A1
CF_PPI_network.lcc.interactors_to_hubs.2 <- CF_PPI_network.lcc.interactors_to_hubs[which(row.names(CF_PPI_network.lcc.interactors_to_hubs)!="CSNK2A1"), which(colnames(CF_PPI_network.lcc.interactors_to_hubs)!="SRC")]

### Hubs to Endpoints
CF_PPI_network.lcc.hubs_to_endpoints <- CF_PPI_network.lcc.dists[hubs, endpoints]

### Hubs to Endpoints
CF_PPI_network.lcc.hubs_to_hubs <- CF_PPI_network.lcc.dists[hubs, hubs]


# Betweenness
CF_PPI_network.lcc.bc <- as.data.frame(betweenness(CF_PPI_network.lcc.igraph))
CF_PPI_network.lcc.bc$HGNC <- rownames(CF_PPI_network.lcc.bc)
colnames(CF_PPI_network.lcc.bc) <- c("Betweeness", "HGNC")
CF_PPI_network.lcc.bc <- CF_PPI_network.lcc.bc[order(CF_PPI_network.lcc.bc$Betweeness, decreasing = T),]
CF_PPI_network.lcc.bc$HGNC <- factor(CF_PPI_network.lcc.bc$HGNC,
                                        levels = CF_PPI_network.lcc.bc$HGNC)
CF_PPI_network.lcc.bc.top20 <- CF_PPI_network.lcc.bc[1:20,]


# Barplot
ggplot(CF_PPI_network.lcc.bc.top20, aes(x=HGNC, y=Betweeness)) +
  geom_bar(stat = "identity")+
  scale_x_discrete(name = "HGNC Symbol")+
  scale_y_continuous(name = "Betweeness")+
  theme(axis.text.x = element_text(face="bold", size=18, angle=90))


# Heatmap

CF_PPI_network.nodes.bc.pathways_connected.df <- CF_PPI_network.lcc.node_type@nodes[which(CF_PPI_network.lcc.node_type@nodes$Symbol %in% as.character(CF_PPI_network.lcc.bc.top20$HGNC)),c("Symbol", diff_pathways)]
CF_PPI_network.nodes.bc.pathways_connected.mat <- as.matrix(CF_PPI_network.nodes.bc.pathways_connected.df[,diff_pathways])
rownames(CF_PPI_network.nodes.bc.pathways_connected.mat) <- CF_PPI_network.nodes.bc.pathways_connected.df$Symbol
CF_PPI_network.nodes.bc.pathways_connected.mat <- CF_PPI_network.nodes.bc.pathways_connected.mat[levels(CF_PPI_network.lcc.bc.top20$HGNC)[1:20],]

bc.heatmap <- pheatmap(t(CF_PPI_network.nodes.bc.pathways_connected.mat),
                         color = c("white", "black"),
                         cluster_rows = F,
                         cluster_cols = F,
                         show_colnames = F,
                         fontsize = 10)


bc.plot <- ggplot(CF_PPI_network.lcc.bc.top20, aes(x=Betwee)) + 
  geom_density()











# End poitns

endpoints.in_CF_PPI_network <- CF_PPI_network.llc.curated.2.node_type@nodes[which(CF_PPI_network.llc.curated.2.node_type@nodes$endpoint_topology),"Symbol"]

endpoints.in_CF_PPI_network.pathways_connected <- lapply(endpoints.in_CF_PPI_network, function(end_point){
  return(pathways_connected_function(end_point, 
                                     CF_PPI_network.llc.curated.2.node_type))
})

endpoints.in_CF_PPI_network.pathways_connected.interactions <- do.call("rbind", 
                                                             endpoints.in_CF_PPI_network.pathways_connected)

endpoints.in_CF_PPI_network.pathways_connected.nodes <- get_network_nodes(network_df = data.frame(endpoints.in_CF_PPI_network.pathways_connected.interactions),
                                                                   source_colname = "genesymbol_source",
                                                                   target_colname = "genesymbol_target")

endpoints.in_CF_PPI_network.pathways_connected.nodes <- merge(endpoints.in_CF_PPI_network.pathways_connected.nodes,
                                                       CF_PPI_network.llc.curated.2.node_type@nodes[,c("Symbol",
                                                                                                       "Node_type", 
                                                                                                       "CFTR_interactor", 
                                                                                                       "Type")],
                                                       by.x = "HGNC",
                                                       by.y = "Symbol",
                                                       all.x=TRUE)

endpoints.in_CF_PPI_network.pathways_connected.nodes[which(is.na(endpoints.in_CF_PPI_network.pathways_connected.nodes$Node_type)),"Node_type"]<- "supernode"

CF_PPI_network.pathways_connected <- rbind(endpoints.in_CF_PPI_network.pathways_connected.interactions,
                                           all_CFTR_interactors.pathways_connected.interactions)

all_CFTR_interactors.pathways_connected.nodes$Count <- NULL
endpoints.in_CF_PPI_network.pathways_connected.nodes$Count <- NULL
CF_PPI_network.pathways_connected.nodes <- unique(rbind(all_CFTR_interactors.pathways_connected.nodes,
                                           endpoints.in_CF_PPI_network.pathways_connected.nodes))


# test.casp <- pathways_connected_function("NFKB1",
#                                     CF_PPI_network.llc.curated.2.node_type)
# apoptosis_nodes <- CF_PPI_network.llc.curated.2.node_type@nodes$Symbol[grep("CASP", CF_PPI_network.llc.curated.2.node_type@nodes$Symbol)]
# 
# apoptosis_interactions <- CF_PPI_network.llc.curated.2.node_type@interactions[which(CF_PPI_network.llc.curated.2.node_type@interactions$genesymbol_source %in% apoptosis_nodes |
#                                                                                       CF_PPI_network.llc.curated.2.node_type@interactions$genesymbol_target %in% apoptosis_nodes),]


# write.table(CF_PPI_network.pathways_connected,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_and_TF_PPathI_interactions_df_2023_03_22.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)
# 
# write.table(CF_PPI_network.pathways_connected.nodes,
#             file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_and_TF_PPathI_nodes_df_2023_03_22.txt",
#             sep = "\t",
#             row.names = F,
#             quote = FALSE)

# write.table(all_CFTR_interactors.interactions.df,
#             "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/CFTR_interactors_interactions_2022_03_14.csv",
#             sep = "\t",
#             col.names = T,
#             row.names = F,
#             quote = F)

all_CFTR_interactors.interactions.nodes <- unique(c(all_CFTR_interactors.interactions.df$genesymbol_source,
                                                    all_CFTR_interactors.interactions.df$genesymbol_target))

all_CFTR_interactors.interactions.nodes.df <- CF_PPI_network.llc.curated.2.node_type@nodes[which(CF_PPI_network.llc.curated.2.node_type@nodes$Symbol %in% all_CFTR_interactors.interactions.nodes),]

# write.table(all_CFTR_interactors.interactions.nodes.df,
#             "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/CFTR_interactors_interactions_nodes_2022_03_14.csv",
#             sep = "\t",
#             col.names = T,
#             row.names = F,
#             quote = F)

# to_test 
protein <- "SYK"
protein <- "STX1A"
network <- CF_PPI_network.llc.curated.2.node_type
pathways_connected_function <- function(protein, network) {
  
  protein.targets_interactions <- network@interactions[which(network@interactions$genesymbol_source==protein),c("effect",
                                                                                                                "arrow",
                                                                                                                "genesymbol_source",
                                                                                                                "genesymbol_target")]
  
  protein.targets <- network@nodes[which(network@nodes$Symbol %in% protein.targets_interactions$genesymbol_target),c("Symbol",
                                                                                                                     diff_pathways,
                                                                                                                     "pathway")]
  
  protein.targets.df <- merge(protein.targets_interactions,
                              protein.targets,
                              by.x = "genesymbol_target",
                              by.y = "Symbol")
  
  protein.targets.df <- protein.targets.df[,c("genesymbol_source",
                                              "genesymbol_target",
                                              "effect",
                                              "arrow",
                                              diff_pathways,
                                              "pathway")]
  # print(protein.targets.df)
  
  protein.sources_interactions <- network@interactions[which(network@interactions$genesymbol_target==protein),c("effect",
                                                                                                                "arrow",
                                                                                                                "genesymbol_source",
                                                                                                                "genesymbol_target")]
  
  protein.sources <- network@nodes[which(network@nodes$Symbol %in% protein.sources_interactions$genesymbol_source),c("Symbol",
                                                                                                                     diff_pathways,
                                                                                                                     "pathway")]
  
  protein.sources.df <- merge(protein.sources_interactions,
                              protein.sources,
                              by.x = "genesymbol_source",
                              by.y = "Symbol")
  
  protein.sources.df <- protein.sources.df[,c("genesymbol_source",
                                              "genesymbol_target",
                                              "effect",
                                              "arrow",
                                              diff_pathways,
                                              "pathway")]
  # print(protein.sources.df)
  
  protein.interactors <- rbind(protein.targets.df,
                               protein.sources.df)
  
  
  # remove interactions with proteins in no pathway
  protein.interactors <- protein.interactors[apply(X = protein.interactors[,diff_pathways],
        MARGIN = 1,
        FUN = function(x){
          return(!all(is.na(x)))
        }),]
  
  # print(protein.interactors)
  
  pathways.connected <- unique(unlist(apply(X = protein.interactors[,diff_pathways], 
                                                                  MARGIN = 1, 
                                                                  FUN = function(x) {
                                                                    return(diff_pathways[x==1])}),
                               use.names = F))
  names(pathways.connected) <- NULL
  pathways.connected.df <- data.frame(CFTR_interactor = protein,
                                      pathway = pathways.connected)
  colnames(pathways.connected.df) <- c("CFTR_interactor",
                                       "pathway")
  
  return(pathways.connected.df)
  
}


all_CFTR_interactors.pathways_connected <- lapply(CFTR_interactor.in_CF_PPI_network$Symbol, function(CFTR.interactor){
  return(pathways_connected_function(CFTR.interactor, 
                                    CF_PPI_network.llc.curated.2.node_type))
})

all_CFTR_interactors.pathways_connected.df <- do.call("rbind", all_CFTR_interactors.pathways_connected)





all_CFTR_interactors.pathways_connected.df$connected <- TRUE
all_CFTR_interactors.pathways_connected.wide <- 
  pivot_wider(data = all_CFTR_interactors.pathways_connected.df,
              names_from = pathway,
              values_from = connected)
all_CFTR_interactors.pathways_connected.wide <- merge(all_CFTR_interactors.pathways_connected.wide,
                                                      CFTR_interactor.in_CF_PPI_network[,c("Symbol",
                                                                                           "Type")],
                                                      by.x = "CFTR_interactor",
                                                      by.y = "Symbol")
all_CFTR_interactors.pathways_connected.wide <- all_CFTR_interactors.pathways_connected.wide[,c("CFTR_interactor",
                                                                                                "Type",
                                                                                                diff_pathways)]

# write.table(all_CFTR_interactors.pathways_connected.wide,
#             "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/CFTR_interactors_pathway_connected_2022_03_14.csv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)

# 2 End points

endpoints.in_CF_PPI_network <- CF_PPI_network.llc.curated.2.node_type@nodes[which(CF_PPI_network.llc.curated.2.node_type@nodes[,"endpoint_to_keep"]),]

endpoints.in_CF_PPI_network <- endpoints.in_CF_PPI_network[,c("Symbol",
                                                              diff_pathways,
                                                              "sum",
                                                              "pathway")]

# write.table(endpoints.in_CF_PPI_network,
#             "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/end_points_pathway_connected_2022_03_15.csv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)