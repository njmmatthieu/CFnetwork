from_complex_to_protein <- function(CyFi_network.interactions,
                                    CyFi_network.nodes,
                                    CFTR_complexes=NULL) {
  
  # Change complex interactions to protein interactions
  rownames(CyFi_network.interactions) <- NULL
  CyFi_network.interactions$Complex.name <- NA 

  # print(CFTR_complexes)
  
  ## REACTANT
  reactant_complex_proteins_nb_all <- list()
  for (i_row in rownames(CyFi_network.interactions)){
  
  # print(i_row)
    if (CyFi_network.interactions[i_row,"REACTANT_Type"]=="Complex" & 
    !(CyFi_network.interactions[i_row,"REACTANT_Name"] %in% CFTR_complexes)){
      reactant_id <- CyFi_network.interactions[i_row, "REACTANT_id"]
      
      # print(reactant_id)
      
      reactant_complex_proteins <- CyFi_network.nodes[which(CyFi_network.nodes$Complex.id==reactant_id),]
      reactant_complex_proteins_nb <- dim(reactant_complex_proteins)[1]
      reactant_complex_proteins_nb_all <- append(reactant_complex_proteins_nb_all,
                                                 reactant_complex_proteins_nb)
      # print(reactant_complex_proteins)
      
      reactant_complex_proteins_interactions <- data.frame(PRODUCT_id=rep(CyFi_network.interactions[i_row,"PRODUCT_id"], reactant_complex_proteins_nb),
                                                           REACTANT_id=reactant_complex_proteins$Id,
                                                           References=rep(CyFi_network.interactions[i_row,"References"], reactant_complex_proteins_nb),
                                                           # Id=rep(CyFi_network.interactions[i_row,"Id"], reactant_complex_proteins_nb),
                                                           Reaction.type=rep(CyFi_network.interactions[i_row,"Reaction.type"], reactant_complex_proteins_nb),
                                                           REACTANT_Name=reactant_complex_proteins$Name,
                                                           REACTANT_Type=reactant_complex_proteins$Type,
                                                           PRODUCT_Name=rep(CyFi_network.interactions[i_row,"PRODUCT_Name"], reactant_complex_proteins_nb),
                                                           PRODUCT_Type=rep(CyFi_network.interactions[i_row,"PRODUCT_Type"], reactant_complex_proteins_nb),
                                                           CFTR_interactions=rep(CyFi_network.interactions[i_row,"CFTR_interactions"], reactant_complex_proteins_nb),
                                                           Complex.name=rep(unique(reactant_complex_proteins$Complex.name),reactant_complex_proteins_nb))
      
      CyFi_network.interactions <- rbind(CyFi_network.interactions,
                                              reactant_complex_proteins_interactions)
      
    }
  }
  if (length(reactant_complex_proteins_nb_all)!=0){
    # print(which(CyFi_network.interactions[,"REACTANT_Type"]=="Complex" & 
    # !(CyFi_network.interactions[i_row,"REACTANT_Name"] %in% CFTR_complexes)))
    CyFi_network.interactions <- CyFi_network.interactions[-which(CyFi_network.interactions[,"REACTANT_Type"]=="Complex" & 
    !(CyFi_network.interactions[,"REACTANT_Name"] %in% CFTR_complexes)),]
  }
  ## PRODUCT
  product_complex_proteins_nb_all <- list()
  for (i_row in rownames(CyFi_network.interactions)){
    
    if (CyFi_network.interactions[i_row,"PRODUCT_Type"]=="Complex" & 
    !(CyFi_network.interactions[i_row,"PRODUCT_Name"] %in% CFTR_complexes)){
      product_id <- CyFi_network.interactions[i_row, "PRODUCT_id"]
      
      # print(product_id)
      
      product_complex_proteins <- CyFi_network.nodes[which(CyFi_network.nodes$Complex.id==product_id),]
      product_complex_proteins_nb <- dim(product_complex_proteins)[1]
      product_complex_proteins_nb_all <- append(product_complex_proteins_nb_all,
                                                product_complex_proteins_nb)
      # print(dim(product_complex_proteins))
      
      product_complex_proteins_interactions <- data.frame(PRODUCT_id=product_complex_proteins$Id,
                                                          REACTANT_id=rep(CyFi_network.interactions[i_row,"REACTANT_id"], product_complex_proteins_nb),
                                                          References=rep(CyFi_network.interactions[i_row,"References"], product_complex_proteins_nb), 
                                                          # Id=rep(CyFi_network.interactions[i_row,"Id"], product_complex_proteins_nb),
                                                          Reaction.type=rep(CyFi_network.interactions[i_row,"Reaction.type"], product_complex_proteins_nb),
                                                          REACTANT_Name=rep(CyFi_network.interactions[i_row,"REACTANT_Name"], product_complex_proteins_nb),
                                                          REACTANT_Type=rep(CyFi_network.interactions[i_row,"REACTANT_Type"], product_complex_proteins_nb),
                                                          PRODUCT_Name=product_complex_proteins$Name,
                                                          PRODUCT_Type=product_complex_proteins$Type,
                                                          CFTR_interactions=rep(CyFi_network.interactions[i_row,"CFTR_interactions"], product_complex_proteins_nb),
                                                          Complex.name=rep(unique(product_complex_proteins$Complex.name),product_complex_proteins_nb))
      
      CyFi_network.interactions <- rbind(CyFi_network.interactions,
                                    product_complex_proteins_interactions)
      
    }
  }
  if (length(product_complex_proteins_nb_all)!=0){
    # print(which(CyFi_network.interactions[,"PRODUCT_Type"]=="Complex" & 
    # !(CyFi_network.interactions[i_row,"PRODUCT_Name"] %in% CFTR_complexes)))
    CyFi_network.interactions <- CyFi_network.interactions[-which(CyFi_network.interactions[,"PRODUCT_Type"]=="Complex" & 
    !(CyFi_network.interactions[,"PRODUCT_Name"] %in% CFTR_complexes)),]
  }
  
  return(CyFi_network.interactions)
}

