library(dorothea)

# for PPI_network class
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/network_utils.R")

# for Endpoints phenotypes
endpoints.pheno.df <- read.table(file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/topological_endpoint_to_phenotype_2023_04_07.txt",
                                 sep = "\t",
                                 header = T,
                                 na.strings = "")
endpoints.pheno.final.df <- endpoints.pheno.df[which(endpoints.pheno.df$Endpoint_bool |
                                                       is.na(endpoints.pheno.df$Endpoint_bool)),]
caspases <- endpoints.pheno.final.df[which(endpoints.pheno.final.df$Endpoint_cat=="Caspases"),"Symbol"]
others <- endpoints.pheno.final.df[which(endpoints.pheno.final.df$Endpoint_cat=="Other"),"Symbol"]
regactin <- endpoints.pheno.final.df[which(endpoints.pheno.final.df$Endpoint_cat=="Regulation of actin cytoskeleton"),"Symbol"]
tf <- endpoints.pheno.final.df[which(endpoints.pheno.final.df$Endpoint_cat=="TF"),"Symbol"]
# # to test
# PPI_network <- CF_PPI_network.lcc.curated.2.weird_endpoint
# source_colname <- "genesymbol_source"
# target_colname <- "genesymbol_target"
# interactors <-  CFTR_interactors
# include_weird_endpoints <- TRUE

# Node type
get_node_type <- function(PPI_network,
                          source_colname="genesymbol_source",
                          target_colname="genesymbol_target",
                          # interactors,
                          include_weird_endpoints=FALSE){
  
  # # Receptors Ligands 
  # receptors_ligands.df <- extract_receptors_ligands(PPI_network)
  # # Receptors
  # receptors.df <- extract_receptors(PPI_network)
  
  PPI_network@nodes$Node_type <- sapply(PPI_network@nodes$Symbol, function(symbol){
    
    # CFTR
    if (symbol=='CFTR'){
      Node_type <- 'CFTR'
    } else if (symbol %in% PPI_network@nodes[which(PPI_network@nodes$receptor_ligand),]$Symbol){
      Node_type <- 'Receptor Ligand'
    } else if (symbol %in% PPI_network@nodes[which(PPI_network@nodes$receptor),]$Symbol){
      Node_type <- 'Receptor'
      # endpoints (TF and Caspase)
    } else if (symbol %in% caspases){
      Node_type <- 'End point - Caspases'
    } else if (symbol %in% others){
    Node_type <- 'End point - Others'
    } else if (symbol %in% regactin){
    Node_type <- 'End point - Regulation of actin cytoskeleton'
    } else if (symbol %in% tf){
    Node_type <- 'End point - TF'
    } else {
      Node_type <- 'Protein'
    }
    
    # # weird endpoint
    # if (include_weird_endpoints) {
    #   if (symbol %in% PPI_network@nodes[which(PPI_network@nodes$weird_endpoint),]$Symbol){
    #     Node_type <- 'Weird endpoint'
    #   }
    # }
    
    # print(Node_type)
    return(Node_type)
    
  })
  
  # PPI_network@nodes$CFTR_interactor <- PPI_network@nodes$Symbol %in% interactors$Symbol 
  # PPI_network@nodes <- merge(PPI_network@nodes,
  #                            interactors,
  #                            by = "Symbol",
  #                            all.x = TRUE)
  
  return(PPI_network)
  
}

## merge logFC of each study

load("/Users/matthieu/ownCloud/Thèse/Systems Biology/Transcriptomic studies/deg_comparison/logFC_combined_2022_06_14.RData")
logFC_combined_df$Symbol <- rownames(logFC_combined_df)

# Node type
get_logFC_studies <- function(PPI_network,
                              source_colname="genesymbol_source",
                              target_colname="genesymbol_target"){
  
  PPI_network@nodes <- merge(PPI_network@nodes,
                             logFC_combined_df, 
                             by = "Symbol",
                             all.x = T)
  
  return(PPI_network)
}