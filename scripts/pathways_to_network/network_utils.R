library(igraph)

setClass("PPI_network", slots=list(interactions="data.frame",
                                         nodes="data.frame"))


preprocess_PPI_network <- function(PPI_network,
                                   source_colname="genesymbol_source",
                                   target_colname="genesymbol_target"){
  
  # Removing IFNA13
  PPI_network@interactions <- PPI_network@interactions[which(PPI_network@interactions$genesymbol_source!="IFNA13" &
                                                               PPI_network@interactions$genesymbol_target!="IFNA13"),]
  PPI_network@nodes <- PPI_network@nodes[which(PPI_network@nodes$Symbol!="IFNA13"),]
  
  # replace DDX58 by RIGI
  PPI_network@nodes$Symbol <- gsub(pattern = "DDX58",
                                   replacement = "RIGI", 
                                   x =PPI_network@nodes$Symbol)
  
  PPI_network@interactions[,source_colname] <- gsub(pattern = "DDX58",
                                                    replacement = "RIGI", 
                                                    x = PPI_network@interactions[,source_colname])
  
  PPI_network@interactions[,target_colname] <- gsub(pattern = "DDX58",
                                                    replacement = "RIGI", 
                                                    x = PPI_network@interactions[,target_colname])
  
  # change effect column
  PPI_network@interactions[which(PPI_network@interactions$effect=="compound"),"effect"] <- paste(
    PPI_network@interactions[which(PPI_network@interactions$effect=="compound"),"effect"],
    PPI_network@interactions[which(PPI_network@interactions$effect=="compound"),"arrow"],
    sep = " - "
  )
  
  return(PPI_network)
  
}


get_network_nodes <- function(network_df,  
                              source_colname="HGNC.A",
                              target_colname="HGNC.B"){
  nodes <- c(network_df[,source_colname], network_df[,target_colname])
  network_nodes <- data.frame(table(nodes))
  colnames(network_nodes) <- c("HGNC", "Count")
  network_nodes$HGNC <- as.character(network_nodes$HGNC)
  return(network_nodes)
}



