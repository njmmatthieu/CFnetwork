library(tibble)

add_kegg_interaction <- function(PPI_network,
                                 Symbol2UniprotID,
                                 source_colname="genesymbol_source",
                                 target_colname="genesymbol_target",
                                 Source,
                                 Effect,
                                 Target,
                                 Pathway_name) {
  
  interactions.df <- PPI_network@interactions
  
  Uniprot_source <- dplyr::filter(Symbol2UniprotID, symbol==Source)[['uniprot_ids']][[1]]
  Uniprot_target <- dplyr::filter(Symbol2UniprotID, symbol==Target)[['uniprot_ids']][[1]]
  Effect_arrow <- effect_arrow.df[which(effect_arrow.df$effect==Effect),"arrow"]
  
  print(Uniprot_source)
  print(Uniprot_target)
  print(Effect_arrow)
  
  new_interaction <- c("PPrel", 
    Effect,
    Effect_arrow,
    Source,
    Uniprot_source,
    Target,
    Uniprot_target,
    Pathway_name,
    # paste(Uniprot_source,
    #       Uniprot_target,
    #       Effect,
    #       sep = "_"),
    TRUE
  )
  
  PPI_network@interactions[nrow(PPI_network@interactions)+1,] <- new_interaction
  print(paste(Source,Effect_arrow, Target,"is added.", sep = " "))

  return(PPI_network)
  
}

add_node <- function(PPI_network,
                     source_colname="genesymbol_source",
                     target_colname="genesymbol_target",
                     Node,
                     Node_UniprotID=NA,
                     Pathway_name) {
  
  new_node <- c(Node,
                Pathway_name)
  
  PPI_network@nodes[nrow(PPI_network@nodes)+1,] <- new_node
  print(paste(Node,"is added.", sep = " "))
  
  return(PPI_network)
  
}

remove_kegg_interaction <- function(PPI_network,
                                 source_colname="genesymbol_source",
                                 target_colname="genesymbol_target",
                                 Source,
                                 Effect,
                                 Target,
                                 Pathway_name) {
  
  
  Effect_arrow <- effect_arrow.df[which(effect_arrow.df$effect==Effect),"arrow"]
  # print(paste("Removing", Source,Effect_arrow, Target, sep = " "))
  
  interaction_id <- which(PPI_network@interactions$genesymbol_source==Source &
                            PPI_network@interactions$genesymbol_target==Target &
                            PPI_network@interactions$effect==Effect)
  
  if (length(interaction_id)==0){
    stop('The interaction does not exist.')
  }
  
  PPI_network@interactions <- PPI_network@interactions[-interaction_id,]
  print(paste(Source,Effect_arrow, Target,"is removed.", sep = " "))
  
  return(PPI_network)
  
}

remove_node <- function(PPI_network,
                                  source_colname="genesymbol_source",
                                  target_colname="genesymbol_target",
                                  Node,
                                  Pathway_name) {
  
  if (!Node %in% PPI_network@nodes$Symbol){
    stop('The protein does not exist.')
  }
  
  PPI_network@nodes <- PPI_network@nodes[which(PPI_network@nodes$Symbol!=Node),]
  interactions_to_remove <- apply(X = PPI_network@interactions[,c(source_colname,target_colname)],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    return(all(x!=Node))
                                  })
  nb_interactions_to_remove <- sum(as.integer(!interactions_to_remove))
  PPI_network@interactions <- PPI_network@interactions[interactions_to_remove,]
  print(paste(Node,"is removed, corresponding to", nb_interactions_to_remove, "interactions.",sep = " "))
  
  return(PPI_network)
  
}

remove_node_uniprotid <- function(PPI_network,
                        UniprotID,
                        Pathway_name) {
  
  interactions_to_remove <- apply(X = PPI_network@interactions[,c("uniprot_source","uniprot_target")],
                                  MARGIN = 1,
                                  FUN = function(x){
                                    return(all(x!=UniprotID))
                                  })
  nb_interactions_to_remove <- sum(as.integer(!interactions_to_remove))
  PPI_network@interactions <- PPI_network@interactions[interactions_to_remove,]
  print(paste(UniprotID,"is removed, corresponding to", nb_interactions_to_remove, "interactions.",sep = " "))
  
  return(PPI_network)
  
}

change_effect <- function(PPI_network,
                                    source_colname="genesymbol_source",
                                    target_colname="genesymbol_target",
                                    Source,
                                    Effect,
                                    Target,
                                    Pathway_name) {
  
  interaction_id <- which(PPI_network@interactions$genesymbol_source==Source &
                            PPI_network@interactions$genesymbol_target==Target)
  if (length(interaction_id)==0){
    stop('The interaction does not exist.')
  }
  previous_effect_arrow <- PPI_network@interactions[interaction_id,"arrow"]
  PPI_network@interactions[interaction_id,"effect"] <- Effect
  
  Effect_arrow <- effect_arrow.df[which(effect_arrow.df$effect==Effect),"arrow"]
  PPI_network@interactions[interaction_id,"arrow"] <- Effect_arrow
  
  print(paste(Source,previous_effect_arrow, Target,"has been changed to",Source,Effect_arrow, Target, sep = " "))
  
  return(PPI_network)
  
}

add_new_nodes_uniprot_ids <- function(Symbol2UniprotID,
                                      corrections.df) {
  corrections.df <- kegg_pathways_corrections
  corrections.df.new_node <- corrections.df[which(corrections.df$action=="add_node"),]
  
  for (i_act in 1:nrow(corrections.df.new_node)){
    # print(i_act)
    # action <- corrections.df[i_act, "action"]
    source <- corrections.df.new_node[i_act, "Source"]
    # print(source)
    source_UniprotID <- corrections.df.new_node[i_act, "Source_UniprotID"]
    # effect <- corrections.df[i_act, "Effect"]
    # target <- corrections.df[i_act, "Target"]
    # pathway_name <- corrections.df[i_act, "pathway_name"]

    if (!source %in% Symbol2UniprotID$symbol){
      print(paste(source, "does not exist in KEGG signaling database.", sep = " "))
      Symbol2UniprotID <- Symbol2UniprotID %>% add_row(symbol = source,
                                                                         uniprot_ids = list(source_UniprotID),
                                                                         gene_group = list(NA))
    }
    
  }
  return(Symbol2UniprotID)
}

# PPI_network <- example_kegg_PPI_network
# Symbol2UniprotID <- kegg_pathways_nodes.carac.corrected
# corrections.df <- kegg_pathways_corrections.example

correct_PPI_network <- function(PPI_network,
                                source_colname="genesymbol_source",
                                target_colname="genesymbol_target",
                                Symbol2UniprotID,
                                corrections.df) {
  
  PPI_network.temp <- PPI_network
  
  for (i_act in 1:nrow(corrections.df)){
    
    action <- corrections.df[i_act, "action"]
    source <- corrections.df[i_act, "Source"]
    source_UniprotID <- corrections.df[i_act, "Source_UniprotID"]
    effect <- corrections.df[i_act, "Effect"]
    target <- corrections.df[i_act, "Target"]
    pathway_name <- corrections.df[i_act, "pathway_name"]
    if (action=="add_interaction"){
      PPI_network.temp <- add_kegg_interaction(PPI_network = PPI_network.temp,
                                               Symbol2UniprotID = Symbol2UniprotID,
                                                            Source = source,
                                                            Effect = effect,
                                                            Target = target,
                                                            Pathway_name = pathway_name)
    } else if (action=="add_node"){
      
      if (!source %in% Symbol2UniprotID$symbol){
        print(paste(source, "does not exist in KEGG signaling database.", sep = " "))
        # Symbol2UniprotID <- Symbol2UniprotID %>% add_row(symbol = source,
        #                                                                    uniprot_ids = list(source_UniprotID),
        #                                                                    gene_group = list(NA))
      }
      
      PPI_network.temp <- add_node(PPI_network = PPI_network.temp,
                                                Node = source,
                                                Pathway_name = pathway_name)
    } else if (action=="remove_interaction"){
      PPI_network.temp <- remove_kegg_interaction(PPI_network = PPI_network.temp,
                                                               Source = source,
                                                               Effect = effect,
                                                               Target = target,
                                                               Pathway_name = pathway_name)
    } else if (action=="remove_node"){
      PPI_network.temp <- remove_node(PPI_network = PPI_network.temp,
                                                   Node = source,
                                                   Pathway_name = pathway_name)
    } else if (action=="remove_node_uniprotid"){
      # print(source_UniprotID)
      PPI_network.temp <- remove_node_uniprotid(PPI_network = PPI_network.temp,
                                                   UniprotID = source_UniprotID,
                                                   Pathway_name = pathway_name)
    } else if (action=="change_effect"){
      # print("lo")
      PPI_network.temp <- change_effect(PPI_network = PPI_network.temp,
                                                     Source = source,
                                                     Effect = effect,
                                                     Target = target,
                                                     Pathway_name = pathway_name)
    }
  }
  
  return(PPI_network.temp)
}