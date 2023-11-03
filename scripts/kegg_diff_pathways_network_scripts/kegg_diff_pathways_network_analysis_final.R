library(igraph)
library(ggplot2)

source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

CF_PPI_network.lcc.node_type.interactions <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_interactions_df_2023_07_10.txt",
                                                        sep = "\t",
                                                        header = T,
                                                        check.names = F)

CF_PPI_network.lcc.node_type.nodes <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathways_with_CFTR_interactors_PPI_lcc_tagged_nodes_df_2023_07_10.txt",
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
upload_GMT_file(gmtfile = "/Users/matthieu/ownCloud/Thèse/Systems Biology/gmt/kegg_from_omnipathR_gsea_2022_09_07.gmt")

CF_PPI_network.nodes.gost.res <- gost(query = CF_PPI_network.lcc.node_type.nodes$Symbol, 
                organism = "gp__uvjl_spkU_zy4", 
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
CF_PPI_network.lcc.node_type.interactions.non_binding <- CF_PPI_network.lcc.node_type@interactions[which(!(CF_PPI_network.lcc.node_type@interactions$effect %in% c("binding/association"))),]
CF_PPI_network.lcc.node_type.interactions.non_binding <- CF_PPI_network.lcc.node_type.interactions.non_binding[,c("genesymbol_source",
                                                                                                                  "genesymbol_target")]
colnames(CF_PPI_network.lcc.node_type.interactions.non_binding) <- c("from", "to")

# ## Dissociation
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
CF_PPI_network.lcc.for_igraph <- rbind(CF_PPI_network.lcc.node_type.interactions.non_binding,
                                       CF_PPI_network.lcc.node_type.interactions.binding.both_directions)

###############
## ENDPOINTS ##
###############

endpoints.final.pathways.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/final_endpoint_to_pathways_2023_04_05_corrected.txt",
                                          sep = "\t",
                                          header = T,
                                          na.strings = "",
                                          check.names = FALSE)
endpoints.final.pathways.df <- endpoints.final.pathways.df[order(endpoints.final.pathways.df$Endpoint_cat),]
endpoints.final.pathways.df$Symbol <- as.character(endpoints.final.pathways.df$Symbol)
endpoints <- as.character(endpoints.final.pathways.df$Symbol)

##################
## Source nodes ##
##################

source_nodes <- c("TRADD",
                  "PRKACA",
                  "SYK",
                  "CSNK2A1",
                  "SRC",
                  "PLCB1",
                  "PLCB3",
                  "EZR")



# Igraph
CF_PPI_network.lcc.igraph <- graph_from_data_frame(CF_PPI_network.lcc.for_igraph, directed=TRUE)

# Degree
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

degree(CF_PPI_network.lcc.igraph, v = CFTR_interactors, mode = "out")



# Betweeness
CF_PPI_network.bc.df <- data.frame(betweenness(CF_PPI_network.lcc.igraph))
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

png(filename="/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/figures/2_4_Analysis_CF_network/betweeness/betweenness_centrality_histogram_2023_07_31.png",
    # units = "cm",
    width=1000,
    height=1000,
    res=100)
CF_PPI_network.bc.hist
dev.off()

CF_PPI_network.bc.df$Symbol <- rownames(CF_PPI_network.bc.df)
rownames(CF_PPI_network.bc.df) <- NULL
colnames(CF_PPI_network.bc.df) <- c("BC.score", "Symbol")

CF_PPI_network.bc.df <- CF_PPI_network.bc.df[order(CF_PPI_network.bc.df$BC.score, decreasing = T),]

write.table(CF_PPI_network.bc.df,
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/betweenness_centrality_CF_network_2023_07_27.csv",
            sep = "\t",
            row.names = F)



## Distances
CF_PPI_network.lcc.dists <- distances(CF_PPI_network.lcc.igraph, mode = "out")

### CFTR interactors to Endpoints
CF_PPI_network.lcc.interactors_to_endpoints <- CF_PPI_network.lcc.dists[source_nodes, endpoints]

endpoints_downstream_to_CFTR_interactors <- data.frame(apply(X = CF_PPI_network.lcc.interactors_to_endpoints,
      MARGIN = 1,
      FUN=function(x){
        return(sum(!is.infinite(x)))
      }))
endpoints_downstream_to_CFTR_interactors$protein <- as.character(rownames(endpoints_downstream_to_CFTR_interactors))
colnames(endpoints_downstream_to_CFTR_interactors) <- c("nb_downstream_endpoints", "protein")
rownames(endpoints_downstream_to_CFTR_interactors) <- NULL
endpoints_downstream_to_CFTR_interactors <- endpoints_downstream_to_CFTR_interactors[order(endpoints_downstream_to_CFTR_interactors$nb_downstream_endpoints,
                                                                                           decreasing = T),]
endpoints_downstream_to_CFTR_interactors$protein <- factor(endpoints_downstream_to_CFTR_interactors$protein,
                         levels=unique(endpoints_downstream_to_CFTR_interactors$protein))

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

barplot_nb_output_nodes <- ggplot(data = endpoints_downstream_to_CFTR_interactors, aes(x = protein, y = nb_downstream_endpoints))+
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

png(filename="/Users/matthieu/ownCloud/Thèse/Systems Biology/Meta-analysis article/figures/2_4_Analysis_CF_network/CFTR_interactors_to_endpoints/nb_downstream_output_nodes_per_protein_2023_10_08.png",
    # units = "cm",
    width=650,
    height=1000,
    res=100)
barplot_nb_output_nodes
dev.off()

## Outputs reached by PRKACA, EZR et CSNK2A1

for (interactor in c("PRKACA", "EZR", "CSNK2A1")){
  # print(interactor)
  # print(CF_PPI_network.lcc.interactors_to_endpoints[interactor,])
  print(colnames(CF_PPI_network.lcc.interactors_to_endpoints)[!is.infinite(CF_PPI_network.lcc.interactors_to_endpoints[interactor,])])
}
View(endpoints.final.pathways.df)

best_candidates <- c("PRKACA", 
                     "EZR", 
                     "CSNK2A1",
                     "PLCB1",
                     "PLCB3",
                     "SYK",
                     "SRC",
                     "TRADD")



# from_candidates_to_selected_toutpus <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/styles/from_4_candidates_to_selected_outputs_2023_06_12.csv",
#                                                   sep = ",",
#                                                   header = T)

CF_network_proteins <- setdiff(CF_PPI_network.lcc.node_type@nodes$Symbol, c(CFTR_interactors, endpoints))

CF_network_proteins.downstream_to_candidates <- CF_PPI_network.lcc.dists[best_candidates, CF_network_proteins]

CF_network_proteins.downstream_to_any_candidate.bool <- apply(X = CF_network_proteins.downstream_to_candidates,
                                                         MARGIN = 2,
                                                         FUN=function(x){
                                                           return(any(!is.infinite(x)))
                                                         })

CF_network_proteins.downstream_to_any_candidate <- colnames(CF_network_proteins.downstream_to_candidates)[CF_network_proteins.downstream_to_any_candidate.bool]


CF_network_proteins.downstream_to_any_candidate.upstream_to_endpoints <- CF_PPI_network.lcc.dists[CF_network_proteins.downstream_to_any_candidate, endpoints]

CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints.bool <- apply(X = CF_network_proteins.downstream_to_any_candidate.upstream_to_endpoints,
                                                              MARGIN = 1,
                                                              FUN=function(x){
                                                                return(all(!is.infinite(x)))
                                                              })

CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints <- rownames(CF_network_proteins.downstream_to_any_candidate.upstream_to_endpoints)[CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints.bool]

CF_PPI_subnetwork.deg.df <- CF_PPI_network.deg.df %>%
  filter(rownames(CF_PPI_network.deg.df) %in% CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints)

CF_PPI_subnetwork.bc.df <- CF_PPI_network.bc.df %>%
  filter(rownames(CF_PPI_network.bc.df) %in% CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints)

CF_PPI_network.bc.df <- data.frame(betweenness(CF_PPI_network.lcc.igraph))
colnames(CF_PPI_network.bc.df) <- "bc.degree"
CF_PPI_network.bc.df$Symbol <- rownames(CF_PPI_network.bc.df) 
  
CF_PPI_network.bc.histo <- ggplot(CF_PPI_network.bc.df, aes(x=bc.degree)) + 
  geom_histogram(binwidth = 500)

CF_PPI_network.bc.df.top <- CF_PPI_network.bc.df[which(CF_PPI_network.bc.df$bc.degree >=4000),]

CF_PPI_network.bc.df.top.upstream_to_endpoints <- CF_PPI_network.lcc.dists[CF_PPI_network.bc.df.top$Symbol, endpoints]

CF_PPI_network.bc.df.top.downstream_to_candidates <- CF_PPI_network.lcc.dists[best_candidates, CF_PPI_network.bc.df.top$Symbol]

CF_PPI_network.bc.df.top <- merge(CF_PPI_network.bc.df.top,
                                  )



CF_subnetwork_proteins.downstream_to_candidates <- CF_PPI_network.lcc.dists[best_candidates,CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints]

CF_PPPI_network_proteins.downstream_to_how_many_candidates <- apply(X = CF_network_proteins.downstream_to_candidates,
                                                              MARGIN = 2,
                                                              FUN=function(x){
                                                                return(sum(!is.infinite(x)))
                                                              })


CF_PPI_network.lcc.node_type.nodes$subnetwork <- CF_PPI_network.lcc.node_type.nodes$Symbol %in% c(best_candidates,
                                                                                                  endpoints,
                                                                                                  CF_network_proteins.downstream_to_any_candidate.upstream_to_all_endpoints.bool)

write.table(CF_PPI_network.lcc.node_type.nodes[,c("Symbol", "subnetwork")],
            file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/diff_kegg_pathway_subnetwork_nodes_2022_07_10.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)

endpoints_downstream_to_int_proteins <- data.frame(apply(X = CF_PPI_network.lcc.proteins_to_endpoints,
                                                             MARGIN = 1,
                                                             FUN=function(x){
                                                               return(sum(!is.infinite(x)))
                                                             }))
colnames(endpoints_downstream_to_int_proteins) <- "distance_to_endpoints"

CF_PPI_network.lcc.CFTR_interactors_to_protein <- CF_PPI_network.lcc.dists[c("TRADD", "SYK", "SRC", "PLCB1"), CF_network_proteins]

CFTR_interactors_to_protein <- data.frame(apply(X = CF_PPI_network.lcc.CFTR_interactors_to_protein,
                                                         MARGIN = 2,
                                                         FUN=function(x){
                                                           return(sum(!is.infinite(x)))
                                                         }))
endpoints_downstream_to_int_proteins$from_CFTR_interactors <- CFTR_interactors_to_protein$apply.X...CF_PPI_network.lcc.CFTR_interactors_to_protein..MARGIN...2..

int_prot <- endpoints_downstream_to_int_proteins[which(endpoints_downstream_to_int_proteins$distance_to_endpoints==35 & endpoints_downstream_to_int_proteins$from_CFTR_interactors==4),]

CF_PPI_network.lcc.CFTR_interactors_to_int_prot <- CF_PPI_network.lcc.dists[c("TRADD", "SYK", "SRC", "PLCB1"), rownames(int_prot)]
CFTR_interactors_to_int_protein <- data.frame(apply(X = CF_PPI_network.lcc.CFTR_interactors_to_int_prot,
                                                MARGIN = 2,
                                                FUN=mean))

CF_PPI_network.lcc.int_proteins_to_endpoints <- CF_PPI_network.lcc.dists[c("TRADD", "SYK", "SRC", "PLCB1"), int_proteins]
CF_PPI_network.lcc.int_proteins_to_endpoints[,"mean"] <- sapply(CF_PPI_network.lcc.int_proteins_to_endpoints, mean)
