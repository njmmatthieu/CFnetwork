library(dorothea)
library(hgnc)
library(tidyverse)

# for PPI_network class, get_network_nodes()
source("scripts/pathways_to_network/pathways_to_network_scripts/network_utils.R")

source("scripts/pathways_to_network/kegg_pathways_manual_curation.R")

# dorotea TF
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

# to test
# PPI_network <- CF_PPI_network.lcc.curated.2
# source_colname <- "genesymbol_source"
# target_colname <- "genesymbol_target"

# Tag dorothea TF and expression interactions
dorothea_tag <-function(PPI_network,
                        source_colname="genesymbol_source",
                        target_colname="genesymbol_target"){
  
  tf_in_kegg_diff_pathways.dorothea <- 
    PPI_network@nodes$Symbol[PPI_network@nodes$Symbol %in% unique(dorothea_regulon_human$tf)]
  
  PPI_network@interactions$expression_dorothea <- 
    (PPI_network@interactions[,source_colname] %in% tf_in_kegg_diff_pathways.dorothea) &  
    !(PPI_network@interactions[,target_colname] %in% tf_in_kegg_diff_pathways.dorothea)
  PPI_network@nodes$dorothea_tf <- 
    PPI_network@nodes$Symbol %in% tf_in_kegg_diff_pathways.dorothea
  
  return(PPI_network)
  
}


# For clarity, remove expression interactions
remove_expression_interactions <- function(PPI_network,
                            source_colname="genesymbol_source",
                            target_colname="genesymbol_target"){

  PPI_network_interactions.not_expression <- 
    PPI_network@interactions[which(!PPI_network@interactions$expression_dorothea & 
                                     !PPI_network@interactions$effect %in% c("expression","repression")),]

  PPI_network_nodes.not_expression <- 
    get_network_nodes(PPI_network_interactions.not_expression,
                      source_colname = "genesymbol_source",
                      target_colname = "genesymbol_target")
  PPI_network_nodes.not_expression <- 
    PPI_network@nodes[which(PPI_network@nodes$Symbol %in% PPI_network_nodes.not_expression$HGNC),]

  PPI_network.without_expression <- new("PPI_network",
                                          interactions=PPI_network_interactions.not_expression,
                                          nodes=PPI_network_nodes.not_expression)

  return(PPI_network.without_expression)

}

# For clarity, remove indirect interactions
# PPI_network <- CF_PPI_network.without_expression
remove_indirect_interactions <- function(PPI_network,
                                           source_colname="genesymbol_source",
                                           target_colname="genesymbol_target"){
  
  PPI_network_interactions.not_indirect <- PPI_network@interactions[which(PPI_network@interactions$effect!="indirect effect"),]
  
  PPI_network_nodes.not_indirect <- get_network_nodes(PPI_network_interactions.not_indirect,
                                                        source_colname = "genesymbol_source",
                                                        target_colname = "genesymbol_target")
  PPI_network_nodes.not_indirect <- PPI_network@nodes[which(PPI_network@nodes$Symbol %in% PPI_network_nodes.not_indirect$HGNC),]
  
  PPI_network.without_indirect <- new("PPI_network",
                                        interactions=PPI_network_interactions.not_indirect,
                                        nodes=PPI_network_nodes.not_indirect)
  
  return(PPI_network.without_indirect)
  
}

####
# Interactions that are activation in one pathway and phophorylation/ubiquitination in another ...
####

same_interactions <- function(PPI_network,
                              source_colname="genesymbol_source",
                              target_colname="genesymbol_target"){
  
  PPI_network@interactions$partial_interaction_id <- apply(X = PPI_network@interactions[,c("uniprot_source",
                                                                                                     "uniprot_target")],
                                                           MARGIN = 1,
                                                           FUN = function(x) {
                                                             # print(x)
                                                             return(paste(x, collapse = "_"))
                                                           })
  duplicated_id <- 
    PPI_network@interactions[duplicated(PPI_network@interactions$partial_interaction_id),
                             "partial_interaction_id" ]
  
  return(PPI_network@interactions[which(PPI_network@interactions$partial_interaction_id %in% duplicated_id),])
}

# to test
# PPI_network <- CF_PPI_network.direct

remove_same_interactions <- function(PPI_network,
                                     source_colname="genesymbol_source",
                                     target_colname="genesymbol_target"){
  
  PPI_network@interactions$partial_interaction_id <- apply(X = PPI_network@interactions[,c("uniprot_source",
                                                                                           "uniprot_target")],
                                                           MARGIN = 1,
                                                           FUN = function(x) {
                                                             # print(x)
                                                             return(paste(x, collapse = "_"))
                                                           })
  duplicated_id <- 
    PPI_network@interactions[duplicated(PPI_network@interactions$partial_interaction_id),
                             "partial_interaction_id" ]
  duplicated_interactions <-  
    PPI_network@interactions[which(PPI_network@interactions$partial_interaction_id %in% duplicated_id),]
  duplicated_interactions.unique <- unique(duplicated_interactions)
  
  duplicated_id.to_remove <- duplicated_id[sapply(duplicated_id, function(id){
    duplicated_effects <-  
      PPI_network@interactions[which(PPI_network@interactions$partial_interaction_id==id),
                               "effect"]
    
    if (any(duplicated_effects=="activation") & any(duplicated_effects=="inhibition")){
      return(FALSE)
    } else {
        return(TRUE)
    }
  })
  ]
  
  id_to_remove <- sapply(duplicated_id.to_remove, function(id){

    interaction_id <- which(PPI_network@interactions$partial_interaction_id==id &
                              PPI_network@interactions$effect!="activation")

    Source <- PPI_network@interactions[interaction_id, source_colname]
    Target <- PPI_network@interactions[interaction_id, target_colname]
    Effect_arrow <- PPI_network@interactions[interaction_id, "arrow"]

    print(paste(Source,Effect_arrow, Target,"is removed.", sep = " "))

    return(interaction_id)
    
  })
  PPI_network@interactions <- PPI_network@interactions[-id_to_remove,]
  PPI_network@interactions$partial_interaction_id <- NULL
  
  return(PPI_network)
  
}

# Remove non source receptors ligands 

receptors_ligands <- c('Receptor ligands',
                       'Chemokine ligands',
                       'Ephrins',
                       'Interferons',
                       'Interleukins',
                       'Tumor necrosis factor superfamily',
                       'Growth hormone family',
                       'VEGF family',
                       'Neurotrophins',
                       'GDNF family ligands',
                       'R-spondin family',
                       'Tachykinin precursors',
                       'Endothelins',
                       'Bone morphogenetic proteins',
                       'Inhibin subunits',
                       'Transforming growth factor beta family')

extract_receptors_ligands <- function(PPI_network,
                            source_colname="genesymbol_source",
                            target_colname="genesymbol_target",
                            genes_groups.df,
                            proteins_to_exclude=c("CFLAR", "RAC1", "CD40LG")){


  nodes_group <- genes_groups.df %>%
  dplyr::filter(symbol %in% PPI_network@nodes$Symbol) %>%
  dplyr::select(c('symbol', 
                  # 'status', 
                  'gene_group'))

  # Receptors Ligands 
  receptors_ligands.df <- nodes_group %>%
    filter(map_lgl(gene_group, ~  any(sapply(receptors_ligands, function(fam) fam %in% .x))))
  receptors_ligands.df <- receptors_ligands.df %>%
    filter(!(symbol %in% proteins_to_exclude))

  return(receptors_ligands.df)

}

# receptors_ligands.df.wide <- receptors_ligands.df %>%
#   +   mutate(gene_group = invoke_map(tibble, gene_group)) %>% 
#   +   unnest()

extract_receptors <- function(PPI_network,
                            source_colname="genesymbol_source",
                            target_colname="genesymbol_target",
                            genes_groups.df){


  nodes_group <- genes_groups.df %>%
  dplyr::filter(symbol %in% PPI_network@nodes$Symbol) %>%
  dplyr::select(c('symbol', 
                  # 'status', 
                  'gene_group'))

  # Receptors
  receptors.df <- nodes_group %>%
    filter(map_lgl(gene_group, ~  any(grepl(pattern = 'receptor', x = .x) &
                                      !(any(sapply(c('Receptor ligands',
                                                     'Estrogen receptors'),
                                                   function(fam) fam %in% .x))))))

  return(receptors.df)

}

tag_prot_cat <- function(PPI_network,
                           source_colname="genesymbol_source",
                           target_colname="genesymbol_target"){
  
  # Import the data set in tidy tabular format
  # NB: Multiple-value columns are kept as list-columns
  (url <- latest_archive_url())
  #> [1] "http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt"
  hgnc_dataset <- import_hgnc_dataset(url)
  # Direct URL to the latest archive in TSV format
  
  # Receptors Ligands 
  receptors_ligands.df <- extract_receptors_ligands(PPI_network, 
                                                    genes_groups.df=hgnc_dataset)
  PPI_network@nodes$receptor_ligand <- PPI_network@nodes$Symbol %in% receptors_ligands.df$symbol
  
  # Receptors
  receptors.df <- extract_receptors(PPI_network, 
                                    genes_groups.df=hgnc_dataset)
  PPI_network@nodes$receptor <- PPI_network@nodes$Symbol %in% receptors.df$symbol
  
  return(PPI_network)
  
}

tag_non_source_receptors_interactions <- function(PPI_network,
                                                  source_colname="genesymbol_source",
                                                  target_colname="genesymbol_target") {

  non_source_receptors <- PPI_network@nodes$Symbol[!PPI_network@nodes$Symbol %in% PPI_network@interactions$genesymbol_source & 
                                                     PPI_network@nodes$receptor]
  
  interactions_to_keep <- apply(X = PPI_network@interactions[,c(source_colname,target_colname)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    return(all(! x %in% non_source_receptors))
                                  })
  
  PPI_network@interactions$to_remove <- TRUE
  PPI_network@interactions[interactions_to_keep,"to_remove"] <- FALSE
  
  nodes_to_keep <- c(PPI_network@interactions[interactions_to_keep,"genesymbol_source"],
                       PPI_network@interactions[interactions_to_keep,"genesymbol_target"])
  
  PPI_network@nodes$to_remove <- TRUE
  PPI_network@nodes[which(PPI_network@nodes$Symbol %in% nodes_to_keep),"to_remove"] <- FALSE
  
  return(PPI_network)
  
}

# # to test
# source_colname="genesymbol_source"
# target_colname="genesymbol_target"
# PPI_network <- CF_PPI_network.curated

remove_non_source_receptors <- function(PPI_network,
                        source_colname="genesymbol_source",
                        target_colname="genesymbol_target") {
  
  non_source_receptors <- PPI_network@nodes$Symbol[!PPI_network@nodes$Symbol %in% PPI_network@interactions$genesymbol_source & 
                                                     PPI_network@nodes$receptor]
  
  interactions_to_keep <- apply(X = PPI_network@interactions[,c(source_colname,target_colname)],
                                MARGIN = 1,
                                FUN = function(x){
                                  return(all(! x %in% non_source_receptors))
                                })
  nb_interactions_to_remove <- sum(as.integer(!interactions_to_keep))
  
  nodes_to_keep <- c(PPI_network@interactions[interactions_to_keep,"genesymbol_source"],
                     PPI_network@interactions[interactions_to_keep,"genesymbol_target"])
  
    
  PPI_network@nodes <- PPI_network@nodes[which(PPI_network@nodes$Symbol %in% nodes_to_keep),]
  PPI_network@interactions <- PPI_network@interactions[interactions_to_keep,]
  
  print(paste(non_source_receptors,"are removed, corresponding to", nb_interactions_to_remove, "interactions.",sep = " "))
  
  return(PPI_network)
  
}

####
# EXTRACT SINGLE CONNECTED COMPONENT
####

# # to test
# source_colname="genesymbol_source"
# target_colname="genesymbol_target"
# PPI_network <- CF_PPI_network.curated.2.w_interactors
# PPI_network <- CF_PPI_network.curated.2

# unique_genes <- unique(c(PPI_network@interactions[,"genesymbol_source"],
#          PPI_network@interactions[,"genesymbol_target"]))

# 
# pathway

# data.frame(pathway_components$membership)

tag_largest_connected_component <- function(PPI_network,
                                            source_colname="genesymbol_source",
                                            target_colname="genesymbol_target"){
  
  pathway_graph <- graph_from_data_frame(PPI_network@interactions[,c("genesymbol_source", 
                                                                     "genesymbol_target")],
                                         directed = TRUE,
                                         vertices = NULL)
  
  pathway_components <- components(graph = pathway_graph,
                                   mode = "weak")
  pathway_components.df <- data.frame(pathway_components$membership)
  pathway_components.df$Symbol <- rownames(pathway_components.df)
  
  largest_component <- which.max(pathway_components$csize)
  pathway_components.df$lcc <- pathway_components.df$pathway_components.membership==largest_component
  
  PPI_network@nodes <- merge(PPI_network@nodes,
                             pathway_components.df,
                             by = "Symbol",
                             all = TRUE)
  
  return(PPI_network)
  
}

# # to test
# source_colname="genesymbol_source"
# target_colname="genesymbol_target"
# PPI_network <- CF_PPI_network.curated

extract_largest_connected_component <- function(PPI_network,
                                             source_colname="genesymbol_source",
                                             target_colname="genesymbol_target"){
  
  pathway_graph <- graph_from_data_frame(PPI_network@interactions[,c("genesymbol_source", 
                                                                 "genesymbol_target")],
                                         directed = TRUE,
                                         vertices = NULL)
  
  pathway_components <- components(graph = pathway_graph,
                                   mode = "weak")
  # pathway_components.df <- data.frame(pathway_components$membership)
  # pathway_components.df$Symbol <- rownames(pathway_components.df)
  # 
  # largest_component <- which.max(pathway_components$csize)
  # pathway_components.df$lcc <- pathway_components.df$pathway_components.membership==largest_component
  # 
  # PPI_network@nodes <- merge(PPI_network@nodes,
  #                            pathway_components.df[,c("Symbol", "lcc")],
  #                            by = "Symbol",
  #                            all = TRUE)
  
  largest_connected_component_protein <- names(pathway_components$membership)[which(pathway_components$membership==which.max(pathway_components$csize))]
  
  largest_connected_component_interactions <- PPI_network@interactions[which(PPI_network@interactions[,source_colname] %in% largest_connected_component_protein &
                                                                            PPI_network@interactions[,target_colname] %in% largest_connected_component_protein),]
  
  
  largest_connected_component_nodes <- PPI_network@nodes %>%
    dplyr::filter(Symbol %in% largest_connected_component_protein)
  # largest_connected_component_nodes <- PPI_network@nodes[which(PPI_network@nodes[,node_colname] %in% largest_connected_component_protein),]
  # print(PPI_network@nodes[,"Symbol"])
  # print(largest_connected_component_protein)
  # print(largest_connected_component_nodes)
  
  largest_connected_component <- new("PPI_network",
                                  interactions=largest_connected_component_interactions,
                                  nodes=largest_connected_component_nodes)
  
  return(largest_connected_component)
  
}


####
# CHECKS
####
# Binding interactions that are in both directions : ProtA --- ProtB and ProtB --- ProtA 
# Not sure that is a real problem
####
binding_interactions <- function(PPI_network,
                                 source_colname="genesymbol_source",
                                 target_colname="genesymbol_target"){
  
  binding_interactions <- PPI_network@interactions[which(PPI_network@interactions$arrow=="---"),]
  both_directions_binding <- which(apply(X = binding_interactions[,c("uniprot_source",
                                                                     "uniprot_target",
                                                                     "effect")],
                                         MARGIN = 1,
                                         FUN = function(x){
                                           return(paste(c(x[2],x[1],x[3]), collapse = "_") %in% binding_interactions$interaction_id)
                                         }))
  both_directions_binding_interactions <- binding_interactions[both_directions_binding,]
  
  return(both_directions_binding_interactions)
}

####
# Weird endpoints
connected.source_or_tf <- function(symbol, 
                                   # gene_to_exclude=NA,
                                   PPI_network) {
  
  # non binding
  PPI_network.non_binding.interactions <- PPI_network@interactions[which(!PPI_network@interactions$arrow %in% c("---", "-+-")),]
  
  interactions.df <- PPI_network@interactions[which(PPI_network@interactions$genesymbol_source==symbol |
                                                      PPI_network@interactions$genesymbol_target==symbol),]
  binding_interactions.df <- interactions.df[which(interactions.df$arrow %in% c("---", "-+-")),]
  
  if(dim(binding_interactions.df)[1]==0){
    return(FALSE)
  } else {
    binding_nodes <- unique(c(binding_interactions.df[,source_colname], 
                              binding_interactions.df[,target_colname]))
    binding_nodes <- setdiff(binding_nodes,
                             symbol)
    binding_nodes.sources <- sapply(binding_nodes, function(binding_symbol){
      return(binding_symbol %in% PPI_network.non_binding.interactions$genesymbol_source)
               # binding_symbol %in% tf)
    })
    return(any(binding_nodes.sources))
  
  }
}