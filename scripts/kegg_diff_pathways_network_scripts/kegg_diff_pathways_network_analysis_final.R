# knitr::opts_knit$set(echo = TRUE, root.dir = normalizePath("../../"))

library(igraph)
library(ggplot2)
library(gprofiler2)

source("scripts/pathways_to_network/network_utils.R")

CF_PPI_network.lcc.node_type.interactions <- 
  read.table(file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df.txt",
             sep = "\t",
             header = T,
             check.names = F)

CF_PPI_network.lcc.node_type.nodes <- 
  read.table(file = "kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df.txt",
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

#############################################
## PATHWAY ENRICHMENT OF THE WHOLE NETWORK ##
#############################################

# "gp__ubAP_9ngg_H0c"
upload_GMT_file(gmtfile = "kegg_pathways/kegg_pathways_from_omnipathR.gmt")

CF_PPI_network.nodes.gost.res <- gost(query = CF_PPI_network.lcc.node_type.nodes$Symbol, 
                organism = "gp__xH8L_h95C_juo", 
                ordered_query = FALSE,
                multi_query = FALSE, 
                significant = FALSE, 
                exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, 
                evcodes = FALSE, 
                user_threshold = 0.05, 
                correction_method = "fdr", # change FDR
                # custom_bg = shared_genes_in_3_studies, 
                # domain_scope = "custom", # change to custom to genes appearing in at least 3 studies
                numeric_ns = "", 
                sources = NULL, 
                as_short_link = FALSE)

############################################################################
## BINDING INTERACTIONS INTO TWO DIRECTED INTERACTIONS IN BOTH DIRECTIONS ##
############################################################################

# Pb of binding interactions
## Non binding
CF_PPI_network.lcc.node_type.interactions.non_binding <- 
  CF_PPI_network.lcc.node_type@interactions[which(!(CF_PPI_network.lcc.node_type@interactions$effect %in% c("binding/association"))),]
CF_PPI_network.lcc.node_type.interactions.non_binding <- 
  CF_PPI_network.lcc.node_type.interactions.non_binding[,c("genesymbol_source",
                                                           "genesymbol_target")]
colnames(CF_PPI_network.lcc.node_type.interactions.non_binding) <- c("from", 
                                                                     "to")

# ## Dissociation
CF_PPI_network.lcc.node_type.interactions.binding <- 
  CF_PPI_network.lcc.node_type@interactions[which(CF_PPI_network.lcc.node_type@interactions$effect %in% c("binding/association")),]
# Both directions for binding interactions
## one direction
CF_PPI_network.lcc.node_type.interactions.binding.one_direction <- 
  CF_PPI_network.lcc.node_type.interactions.binding[,c("genesymbol_source",
                                                       "genesymbol_target")]
colnames(CF_PPI_network.lcc.node_type.interactions.binding.one_direction) <- c("from", 
                                                                               "to")
## other direction
CF_PPI_network.lcc.node_type.interactions.binding.other_direction <- 
  CF_PPI_network.lcc.node_type.interactions.binding[,c("genesymbol_target", 
                                                       "genesymbol_source")]
colnames(CF_PPI_network.lcc.node_type.interactions.binding.other_direction) <- c("from", 
                                                                                 "to")
## both directions
CF_PPI_network.lcc.node_type.interactions.binding.both_directions <- rbind(CF_PPI_network.lcc.node_type.interactions.binding.one_direction,
                                                                           CF_PPI_network.lcc.node_type.interactions.binding.other_direction)
CF_PPI_network.lcc.for_igraph <- rbind(CF_PPI_network.lcc.node_type.interactions.non_binding,
                                       CF_PPI_network.lcc.node_type.interactions.binding.both_directions)

########################################
## SINK NODES AS DEFINED IN THE ARTICLE##
########################################

sink_nodes.final.pathways.df <- read.table(file = "sink_nodes/CFnetwork_sink_nodes_to_pathways.txt",
                                          sep = "\t",
                                          header = T,
                                          na.strings = "",
                                          check.names = FALSE)
sink_nodes.final.pathways.df <- sink_nodes.final.pathways.df[order(sink_nodes.final.pathways.df$Endpoint_cat),]
sink_nodes.final.pathways.df$Symbol <- as.character(sink_nodes.final.pathways.df$Symbol)
sink_nodes <- as.character(sink_nodes.final.pathways.df$Symbol)

##################
## SOURCE NODES ##
##################

source_nodes <- c("TRADD",
                  "PRKACA",
                  "SYK",
                  "CSNK2A1",
                  "SRC",
                  "PLCB1",
                  "PLCB3",
                  "EZR")


#################
## 1. ANALYSIS ##
#################

# Igraph
CF_PPI_network.lcc.igraph <- graph_from_data_frame(CF_PPI_network.lcc.for_igraph, 
                                                   directed=TRUE)

#################
## 1.1 DEGREE ###
#################

CF_PPI_network.deg.df <- data.frame(degree(CF_PPI_network.lcc.igraph, mode = "out"))
CF_PPI_network.deg.hist <-ggplot(CF_PPI_network.deg.df, 
                                 aes(x=degree.CF_PPI_network.lcc.igraph..mode....out..)) + 
  geom_histogram(color="#EA6B66", fill="#F19C99", binwidth = 1)+
  xlab("Degree Centrality")+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=24),
        axis.text.x = element_text(size=24, 
                                   # angle = 45, 
                                   vjust = 4,
                                   hjust=0),
        axis.text.y = element_text(size=24,
                                   hjust=1),
        axis.ticks = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# degree(CF_PPI_network.lcc.igraph, v = CFTR_interactors, mode = "out")

#################################
## 1.1 BETWEENNESS CENTRALITY ###
#################################

CF_PPI_network.bc.df <- data.frame(betweenness(CF_PPI_network.lcc.igraph))
CF_PPI_network.bc.df$Symbol <- rownames(CF_PPI_network.bc.df)
rownames(CF_PPI_network.bc.df) <- NULL
colnames(CF_PPI_network.bc.df) <- c("BC.score", "Symbol")
CF_PPI_network.bc.df <- CF_PPI_network.bc.df[order(CF_PPI_network.bc.df$BC.score, 
                                                   decreasing = T),]

CF_PPI_network.bc.hist <-ggplot(CF_PPI_network.bc.df, 
                                 aes(x=BC.score)) + 
  geom_histogram(fill="#8DA0CB", color="#8DA0CB", binwidth = 500)+
  xlab("Betweenness Centrality")+
  # scale_x_continuous(+
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size=24),
        axis.text.x = element_text(size=24,
                                   angle = 45,
                                   vjust = 0.5,
                                   hjust=0),
        axis.text.y = element_text(size=24,
                                   hjust=1.5),
        axis.ticks = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  scale_x_continuous(breaks = c(0, 2000, 4000, 6000, 8000))

####################
## 1.1 DISTANCES ###
####################

CF_PPI_network.lcc.dists <- distances(CF_PPI_network.lcc.igraph, mode = "out")

### CFTR interactors to sink nodes
CF_PPI_network.lcc.interactors_to_sink_nodes <- CF_PPI_network.lcc.dists[source_nodes, sink_nodes]

sink_nodes_downstream_to_CFTR_interactors <- data.frame(apply(X = CF_PPI_network.lcc.interactors_to_sink_nodes,
      MARGIN = 1,
      FUN=function(x){
        return(sum(!is.infinite(x)))
      }))
sink_nodes_downstream_to_CFTR_interactors$protein <- as.character(rownames(sink_nodes_downstream_to_CFTR_interactors))
colnames(sink_nodes_downstream_to_CFTR_interactors) <- c("nb_downstream_sink_nodes", "protein")
rownames(sink_nodes_downstream_to_CFTR_interactors) <- NULL
sink_nodes_downstream_to_CFTR_interactors <- sink_nodes_downstream_to_CFTR_interactors[order(sink_nodes_downstream_to_CFTR_interactors$nb_downstream_sink_nodes,
                                                                                           decreasing = T),]
sink_nodes_downstream_to_CFTR_interactors$protein <- factor(sink_nodes_downstream_to_CFTR_interactors$protein,
                         levels=unique(sink_nodes_downstream_to_CFTR_interactors$protein))

which(names(V(CF_PPI_network.lcc.igraph))=="PYCARD")
# 322

which(names(V(CF_PPI_network.lcc.igraph))=="NFKB1")
# 322

shortest_paths_from_prot <- function(source, target) {
  
  source_int <- which(names(V(CF_PPI_network.lcc.igraph))==source)
  target_int <- which(names(V(CF_PPI_network.lcc.igraph))==target)
  
  print(all_shortest_paths(CF_PPI_network.lcc.igraph,
                           from=V(CF_PPI_network.lcc.igraph)[source_int],
                           to=V(CF_PPI_network.lcc.igraph)[target_int],
                           mode="out")$res)
}


## Figure

barplot_nb_output_nodes <- ggplot(data = sink_nodes_downstream_to_CFTR_interactors, aes(x = protein, y = nb_downstream_sink_nodes))+
  geom_bar(stat="identity", fill="#CCCCCC") + 
  # ggtitle("Number of downstream \n output nodes \n per protein") +
  ylab("Nb output nodes")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=24, 
                                   angle = 45, 
                                   vjust = 1,
                                   hjust=1),
        axis.text.y = element_text(size=18),
        axis.ticks = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())