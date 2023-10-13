# load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_nodes_df_2023_03_30.RData")
# load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_interactions_df_2023_03_30.RData")

# library(biomaRt)

# get_node_type.CFTR_interactors <- function(PPI_network,
#                                            source_colname="genesymbol_source",
#                                            target_colname="genesymbol_target",
#                                            interactors){
#   
#   PPI_network@nodes$CFTR_interactor <- PPI_network@nodes$Symbol %in% interactors$Symbol 
#   PPI_network@nodes <- merge(PPI_network@nodes,
#                              interactors,
#                              by = "Symbol",
#                              all.x = TRUE)
#   # 
#   # # Node type 2
#   # PPI_network@nodes$Node_type.2 <- PPI_network@nodes$Node_type
#   # PPI_network@nodes[which(PPI_network@nodes$Symbol %in% interactors$Symbol),"Node_type.2"]<- 
#   #   PPI_network@nodes[which(PPI_network@nodes$Symbol %in% interactors$Symbol),"Type"]
#   
#   return(PPI_network)
#   
# }

add_node.CFTR_interactor <- function(PPI_network,
                                     source_colname="genesymbol_source",
                                     target_colname="genesymbol_target",
                                     Node,
                                     Node_UniprotID=NA) {
  
  # PPI_network@nodes$status <- NA
  PPI_network@interactions$reference <- 'KEGG'
  
  PPI_network@nodes[nrow(PPI_network@nodes)+1,"Symbol"] <- Node
  print(paste(Node,"is added.", sep = " "))
  
  return(PPI_network)
  
}

add_kegg_interaction.CFTR_interactor <- function(PPI_network,
                                                 Symbol2UniprotID = kegg_pathways_nodes.carac.corrected,
                                                 source_colname="genesymbol_source",
                                                 target_colname="genesymbol_target",
                                                 Source,
                                                 Effect,
                                                 Target,
                                                 Status) {
  
  # # to test 
  # PPI_network <- PPI_network.temp
  # i_row <- 1 
  # Source <- interactions.str_connected.to_add.df[i_row, "Source"]
  # Effect <- interactions.str_connected.to_add.df[i_row, "Effect"]
  # Target <- interactions.str_connected.to_add.df[i_row, "Target"]
  # Status <- interactions.str_connected.to_add.df[i_row, "status"]
  
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  
  # Ajouter une colonne pour dire que ca vient de KEGG et l'autre que ca vient de la littérature
  if (!all(c("status", "reference") %in% colnames(PPI_network@interactions))){
    PPI_network@interactions$status <- NA
    PPI_network@interactions$reference <- 'KEGG' 
  }
  
  print("coucou")
  
  if (Source %in% Symbol2UniprotID[,"symbol"]){
    print(Source)
    Uniprot_source <- dplyr::filter(Symbol2UniprotID, symbol==Source)[['uniprot_ids']][[1]]
  } else {
  #   Source.id <- getBM(filters= "hgnc_symbol", 
  #         attributes= c("ensembl_gene_id","uniprotswissprot", "hgnc_symbol"),
  #         values=Source,
  #         mart= mart)
  #   Uniprot_source <- setdiff(unique(Source.id$uniprotswissprot), "")
    Uniprot_source <- dplyr::filter(all_CFTR_interactions.PPI.connected.nodes.df, HGNC==Source)[['UniProtID']][[1]]
  }
  if (Target %in% Symbol2UniprotID[,"symbol"]){
    print(Target)
    Uniprot_target <- dplyr::filter(Symbol2UniprotID, symbol==Target)[['uniprot_ids']][[1]]
  } else {
  #   Target.id <- getBM(filters= "hgnc_symbol", 
  #                      attributes= c("ensembl_gene_id","uniprotswissprot", "hgnc_symbol"),
  #                      values=Target,
  #                      mart= mart)
  #   Uniprot_target <- setdiff(unique(Target.id$uniprotswissprot), "")
    Uniprot_target <- dplyr::filter(all_CFTR_interactions.PPI.connected.nodes.df, HGNC==Target)[['UniProtID']][[1]]
  }
  
  
  Effect_arrow <- effect_arrow.df[which(effect_arrow.df$effect==Effect),"arrow"]
  
  print(Uniprot_source)
  print(Uniprot_target)
  print(Effect_arrow)
  
  interaction_id <- paste(Uniprot_source,
                          Uniprot_target,
                          Effect,
                          sep = "_")
  
  new_interaction <- c(interaction_id,
                       Effect,
                       Effect_arrow,
                       Source,
                       Uniprot_source,
                       Target,
                       Uniprot_target,
                       'CyFi-MAP',
                       Status)
  
  PPI_network@interactions[nrow(PPI_network@interactions)+1,c('interaction_id',
                                                              'effect',
                                                              'arrow',
                                                              'genesymbol_source',
                                                              'uniprot_source',
                                                              'genesymbol_target',
                                                              'uniprot_target',
                                                              'reference',
                                                              'status')] <- new_interaction
  print(paste(Source,Effect_arrow, Target,"is added.", sep = " "))
  
  return(PPI_network)
  
}

# # to test
# PPI_network <- CF_PPI_network

extend_to_CFTR_interactors <- function(PPI_network,
                                       include_CFTR=TRUE){
  
  PPI_network.temp <- PPI_network
  
  # 1 - Adding interactors
  
  # interactors inside the network
  all_CFTR_interactors <- all_CFTR_interactions.PPI.connected.nodes.df[which(all_CFTR_interactions.PPI.connected.nodes.df$CFTR_interactor),"HGNC"]
  all_CFTR_interactors.in_network <- all_CFTR_interactors[all_CFTR_interactors %in% PPI_network.temp@nodes$Symbol]
  
  # interactors connetced to the
  interactors.str_connected.to_add.bool <- sapply(all_CFTR_interactors, function(symbol){
    
    # # to test
    # symbol <- "PRKACA"
    # symbol <- "CSNK2A1"
    # symbol <- "SYK"
    # symbol <- "TRADD"
    # symbol <- "HSP90AA1"
    
    interactors.interactions <- all_CFTR_interactions.PPI.connected.min[which(all_CFTR_interactions.PPI.connected.min$Source==symbol |
                                                                                all_CFTR_interactions.PPI.connected.min$Target==symbol),]
    interactors.str_connected <- setdiff(c(interactors.interactions$Source,
                                           interactors.interactions$Target),
                                         c(symbol, "CFTR"))
    
    interactors.in_network <- interactors.str_connected[interactors.str_connected %in% setdiff(PPI_network.temp@nodes$Symbol,all_CFTR_interactors.in_network)]
    # interactors.in_network <- interactors.str_connected[interactors.str_connected %in% PPI_network.temp@nodes$Symbol]

    
    
    return(length(interactors.in_network)>0)
  })
  interactors.str_connected.to_add <- names(interactors.str_connected.to_add.bool)[interactors.str_connected.to_add.bool]
  interactors.connected.to_add <- c(interactors.str_connected.to_add, all_CFTR_interactors.in_network)
  
  if (length(interactors.connected.to_add)==0) {
    stop("No interactor to connect to this pathway")
  }
    
  for (interactor in interactors.connected.to_add){
    
    
    if (!interactor %in% PPI_network.temp@nodes$Symbol){
    
      PPI_network.temp <- add_node.CFTR_interactor(PPI_network = PPI_network.temp,
                                                 Node = interactor)
    }
  }
  
  if (include_CFTR){
    PPI_network.temp <- add_node.CFTR_interactor(PPI_network = PPI_network.temp,
                                                 Node = 'CFTR')
  }
  
  
  # 2 - Adding interactions
  # interactions_to_add.nodes <- unique(c(all_CFTR_interactors, interactors.str_connected.to_add))
  interactions.str_connected.to_add.list <- lapply(interactors.connected.to_add, function(symbol){
    
    print(symbol)
    # # to test 
    # symbol <- "TRADD"
    
    interactors.interactions <- all_CFTR_interactions.PPI.connected.min[which(all_CFTR_interactions.PPI.connected.min$Source==symbol |
                                                                                all_CFTR_interactions.PPI.connected.min$Target==symbol),]
    interactors.str_connected <- setdiff(c(interactors.interactions$Source,
                                           interactors.interactions$Target),
                                         # c(symbol, "CFTR"))
                                         c("CFTR"))
    
    
    # interactors.in_network <- interactors.str_connected[interactors.str_connected%in% PPI_network@nodes$Symbol]
    interactors.in_network <- interactors.str_connected[interactors.str_connected%in%  setdiff(PPI_network.temp@nodes$Symbol,all_CFTR_interactors.in_network)]
    
    if (length(interactors.in_network)>0 | symbol %in% all_CFTR_interactors.in_network){
      interactions.str_connected.to_add <- all_CFTR_interactions.PPI.connected.min[which(apply(X = all_CFTR_interactions.PPI.connected.min[,c("Source", "Target")],
                                                                                               MARGIN = 1,
                                                                                               FUN = function(x){return(all(x %in% c(interactors.in_network, symbol, "CFTR")))})),]
      
      return(interactions.str_connected.to_add)
    } else {
      return(NA)
    }
  })
  interactions.str_connected.to_add.list <- interactions.str_connected.to_add.list[!is.na(interactions.str_connected.to_add.list)]
  interactions.str_connected.to_add.df <- do.call(rbind, interactions.str_connected.to_add.list)
  
  if (!include_CFTR){
    interactions.str_connected.to_add.df <- interactions.str_connected.to_add.df[!interactions.str_connected.to_add.df$CFTR_interactions,]
  }
  print('status' %in% colnames(PPI_network.temp@nodes))
  
  interactions.str_connected.to_add.df <- unique(interactions.str_connected.to_add.df)
  
  PPI_network.temp@interactions$CyFi_MAP_interactions <- FALSE
  
  for (i_row in 1:nrow(interactions.str_connected.to_add.df)){
    
    # i_row <- 1 

    source <- interactions.str_connected.to_add.df[i_row, "Source"]
    effect <- interactions.str_connected.to_add.df[i_row, "Effect"]
    target <- interactions.str_connected.to_add.df[i_row, "Target"]
    status <- interactions.str_connected.to_add.df[i_row, "status"]
    
    print(source)
    print(target)
    
    PPI_network.temp <- add_kegg_interaction.CFTR_interactor(PPI_network = PPI_network.temp,
                                                             Source = source,
                                                             Effect = effect,
                                                             Target = target,
                                                             Status = status)
  }
  
  PPI_network.temp@interactions <- unique(PPI_network.temp@interactions)
  PPI_network.temp@interactions[which(is.na(PPI_network.temp@interactions$CyFi_MAP_interactions)),"CyFi_MAP_interactions"] <- TRUE
  
  # add tag CFTR interactor
  # PPI_network.temp@nodes$CFTR_interactor <- PPI_network.temp@nodes$Symbol %in% all_CFTR_interactors
  PPI_network.temp@nodes <- merge(PPI_network.temp@nodes,
                                  all_CFTR_interactions.PPI.connected.nodes.df[,c("HGNC", 
                                                                                  "CFTR_interactor")], 
                                                                                  # "status")],
                                  all.x = TRUE,
                                  by.x="Symbol",
                                  by.y = "HGNC")
  
  return(PPI_network.temp)
  
}
