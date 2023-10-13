library(dplyr)
library(stringr)

library(biomaRt)

# for from_complex_to_protein
source("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/pathways_to_network_scripts/CyFi_MAP_helper.R")

#############
# WILD TYPE #
#############

## Interactome
WT_interactome <- read.delim("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/interactome/CyFi-MAP/CFTR_cp17-elementExport_2022_10_26_MN_2.txt",
                             sep = "\t",
                             header = T)
## Interactions
WT_interactions <- read.delim("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/interactome/CyFi-MAP/CFTR_cp17-networkExport_2022_10_26.txt",
                              sep = "\t",
                              header = T)

## 1. PREPROCESSING

### Remove empty symbols
WT_interactome[which(WT_interactome$Symbol==""),"Symbol"] <- NA

### WARNING: PB with TRADD
WT_interactions[which(WT_interactions$References=="PUBMED:31231217,PUBMED:27960153"),"Reaction.type"] <- "Negative influence"

### Removing na columns
WT_interactions[colnames(WT_interactions)[sapply(WT_interactions, function(x) return(all(is.na(x))))]]<- NULL
WT_interactions[,c("Id", "Symbol", "Abbreviation", "Synonyms", "Reaction.external.id","Map.id","Map.name")]<- NULL

### replace 'Positive influence' and 'Negative influence' into 'activation' and 'inhibition'
WT_interactions$Reaction.type <- as.factor(WT_interactions$Reaction.type)
levels(WT_interactions$Reaction.type) <- list(activation="Positive influence", inhibition="Negative influence")

### decomposing 'Elements'
WT_interactions$REACTANT_id <- sapply(WT_interactions$Elements, function(reaction){
  return(str_match(string = reaction, pattern = "^REACTANT:([0-9]{7})")[2])
}) 
WT_interactions$PRODUCT_id <- sapply(WT_interactions$Elements, function(reaction){
  return(str_match(string = reaction, pattern = "^REACTANT:[0-9]{7},PRODUCT:([0-9]{7})$")[2])
})
WT_interactions$Elements <- NULL

### Adding Gene Names to interactions dataframe
WT_interactions <- merge(WT_interactions,
                         WT_interactome[,c("Id","Name", "Type")],
                         by.x = "REACTANT_id",
                         by.y = "Id",
                         all.x = TRUE)
WT_interactions <- merge(WT_interactions,
                         WT_interactome[,c("Id","Name", "Type")],
                         by.x = "PRODUCT_id",
                         by.y = "Id",
                         all.x = TRUE)
colnames(WT_interactions)[(ncol(WT_interactions)-3):ncol(WT_interactions)] <- c("REACTANT_Name",
                                                                                "REACTANT_Type",
                                                                                "PRODUCT_Name",
                                                                                "PRODUCT_Type")
## 2. KEEPING INTERACTIONS AT THE PM

WT_interactome.PM.symbol <- WT_interactome[which(WT_interactome$Plasma.Membrane), "Name"]
WT_interactome.not_PM.symbol <- WT_interactome[which(!WT_interactome$Plasma.Membrane), "Name"]

WT_interactome.not_PM.symbol.unique <- unique(setdiff(WT_interactome.not_PM.symbol,
                                                      "CFTR"))
WT_interactions$not_PM <- apply(X = WT_interactions[,c("REACTANT_Name",
                                                                    "PRODUCT_Name")],
                                             MARGIN = 1,
                                             FUN = function(x){
                                               return(any(x %in% WT_interactome.not_PM.symbol.unique))
                                             })
WT_interactions <- WT_interactions[which(!WT_interactions$not_PM),]
WT_interactions[,c("not_PM")]<- NULL

## 3. TAG CFTR INTERACTIONS

WT_interactions$CFTR_interactions <- apply(X = WT_interactions[,c("REACTANT_Name",
                                                                  "PRODUCT_Name")],
                                           MARGIN = 1,
                                           FUN = function(x){
                                             return(any(x=="CFTR") & all(!is.na(x)))
                                           })

## 4. FROM COMPLEX TO PROTEINS

WT_interactions.PPI <- from_complex_to_protein(CyFi_network.interactions = WT_interactions,
                                               CyFi_network.nodes = WT_interactome)

## 5. REMOVE CFTR --- CFTR INTERACTIONS

WT_interactions.PPI <- WT_interactions.PPI[apply(WT_interactions.PPI[,c("REACTANT_Name", "PRODUCT_Name")], 
                                                 MARGIN = 1, 
                                                 FUN = function(x){return(any(x!='CFTR'))}),]

## 6. REMOVE INTERACTIONS WITH 'SIMPLE MOLECULE' AND 'DEGRADED'

WT_interactions.PPI <- WT_interactions.PPI[apply(WT_interactions.PPI[,c("REACTANT_Type", "PRODUCT_Type")], 
                                                 MARGIN = 1, 
                                                 FUN = function(x){return(all(!(x %in% c('Simple molecule', 'Degraded'))))}),]


## 7. KEEP INTERACTIONS WITH 1 INTERMEDIATE TO CFTR

### Direct interactors
WT_CFTR_interactions.PPI <- WT_interactions.PPI[which(WT_interactions.PPI$CFTR_interactions),]
WT_CFTR_interactors <- setdiff(unique(c(WT_CFTR_interactions.PPI$REACTANT_Name, 
                                        WT_CFTR_interactions.PPI$PRODUCT_Name)),
                               "CFTR")
WT_CFTR_interactions.PPI.connected <- WT_interactions.PPI[apply(X=WT_interactions.PPI[,c("REACTANT_Name", "PRODUCT_Name")], 
                                                                MARGIN = 1, 
                                                                FUN = function(x){
                                                                  return(any(x %in% WT_CFTR_interactors))
                                                                }),]

########
# F508 #
########

## Interactome

F508del_interactome <- read.delim("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/interactome/CyFi-MAP/F508del_cp21-elementExport_MN.txt",
                                  sep = "\t",
                                  header = T)

## Interactions

F508del_interactions <- read.delim("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/interactome/CyFi-MAP/F508del_cp21-networkExport.txt",
                                   sep = "\t",
                                   header = T)

# it will be used in every part of the script
F508del_CFTR <- c("F508del",
                  "rF508del",
                  "'Chaperone trap'",
                  "F508del-DNAJB1-HSPA4-HSP90")

## 1. PREPROCESSING

### Remove empty symbols
F508del_interactome[which(F508del_interactome$Symbol==""),"Symbol"] <- NA

### Removing na columns
F508del_interactions[colnames(F508del_interactions)[sapply(F508del_interactions, function(x) return(all(is.na(x))))]]<- NULL
F508del_interactions[,c("Id", "Reaction.external.id", "Map.id","Map.name")]<- NULL

### replace 'Positive influence' and 'Negative influence' into 'activation' and 'inhibition'
F508del_interactions$Reaction.type <- as.factor(F508del_interactions$Reaction.type)
levels(F508del_interactions$Reaction.type) <- list(activation="Positive influence", inhibition="Negative influence")

### decomposing 'Elements'
F508del_interactions$REACTANT_id <- sapply(F508del_interactions$Elements, function(reaction){
  return(str_match(string = reaction, pattern = "^REACTANT:([0-9]{7})")[2])
}) 
F508del_interactions$PRODUCT_id <- sapply(F508del_interactions$Elements, function(reaction){
  return(str_match(string = reaction, pattern = "^REACTANT:[0-9]{7},PRODUCT:([0-9]{7})$")[2])
})
F508del_interactions$Elements <- NULL

# Adding Gene Names to interactions dataframe
F508del_interactions <- merge(F508del_interactions,
                              F508del_interactome[,c("Id","Name", "Type")],
                              # F508del_interactome[,c("Id","Name", "Type","Plasma.Membrane")], # if "Plama.Membrane" column is manually added
                              by.x = "REACTANT_id",
                              by.y = "Id",
                              all.x = TRUE)

F508del_interactions <- merge(F508del_interactions,
                              F508del_interactome[,c("Id","Name", "Type")],
                              # F508del_interactome[,c("Id","Name", "Type","Plasma.Membrane")], # if "Plama.Membrane" column is manually added
                              by.x = "PRODUCT_id",
                              by.y = "Id",
                              all.x = TRUE)
colnames(F508del_interactions)[(ncol(F508del_interactions)-3):ncol(F508del_interactions)] <- c("REACTANT_Name",
                                                                                               "REACTANT_Type",
                                                                                               "PRODUCT_Name",
                                                                                               "PRODUCT_Type")

## 2. KEEPING INTERACTIONS AT THE PM

F508del_interactome.PM.symbol <- F508del_interactome[which(F508del_interactome$Plasma.Membrane), "Name"]
F508del_interactome.not_PM.symbol <- F508del_interactome[which(!F508del_interactome$Plasma.Membrane), "Name"]
F508del_interactome.not_PM.symbol.unique <- unique(setdiff(F508del_interactome.not_PM.symbol,
                                                           F508del_CFTR))

F508del_interactions$not_PM <- apply(X = F508del_interactions[,c("REACTANT_Name",
                                                                 "PRODUCT_Name")],
                                                  MARGIN = 1,
                                                  FUN = function(x){
                                                    return(any(x %in% F508del_interactome.not_PM.symbol.unique))
                                                  })

F508del_interactions <- F508del_interactions[which(!F508del_interactions$not_PM),]
F508del_interactions[,c("not_PM")]<- NULL
F508del_interactome$not_PM <- NULL

## 3. TAG CFTR INTERACTIONS

F508del_interactions$CFTR_interactions <- apply(X = F508del_interactions[,c("REACTANT_Name",
                                                                            "PRODUCT_Name")],
                                                MARGIN = 1,
                                                FUN = function(x){
                                                  return(any(x %in% F508del_CFTR) & 
                                                           all(!is.na(x)))})

## 4. FROM COMPLEX TO PROTEINS

F508del_complexes <- c("'Chaperone trap'",
                       "F508del-DNAJB1-HSPA4-HSP90")

F508del_interactions.PPI <- from_complex_to_protein(CyFi_network.interactions = F508del_interactions,
                                                    CyFi_network.nodes = F508del_interactome,
                                                    CFTR_complexes = F508del_complexes)

## 5. REMOVE CFTR --- CFTR INTERACTIONS

F508del_interactions.PPI <- F508del_interactions.PPI[apply(F508del_interactions.PPI[,c("REACTANT_Name", "PRODUCT_Name")], 
                                                           MARGIN = 1, 
                                                           FUN = function(x){return(any(!(x %in% F508del_CFTR)))}),]

## 6. REMOVE INTERACTIONS WITH 'SIMPLE MOLECULE' AND 'DEGRADED'

F508del_interactions.PPI <- F508del_interactions.PPI[apply(F508del_interactions.PPI[,c("REACTANT_Type", "PRODUCT_Type")], 
                                                           MARGIN = 1, 
                                                           FUN = function(x){return(all(!(x %in% c('Simple molecule', 'Degraded'))))}),]

## 7. KEEP INTERACTIONS WITH 1 INTERMEDIATE TO CFTR

# Direct interactors 
F508del_CFTR_interactions.PPI <- F508del_interactions.PPI[which(F508del_interactions.PPI$CFTR_interactions),]
F508del_CFTR_interactors <- setdiff(unique(c(F508del_CFTR_interactions.PPI$REACTANT_Name, F508del_CFTR_interactions.PPI$PRODUCT_Name)),
                                    F508del_CFTR)
F508del_CFTR_interactions.PPI.connected <- F508del_interactions.PPI[apply(X=F508del_interactions.PPI[,c("REACTANT_Name", "PRODUCT_Name")], 
                                                                          MARGIN = 1, 
                                                                          FUN = function(x){
                                                                            return(any(x %in% F508del_CFTR_interactors))
                                                                          }),]

######################
# MERGE INTERACTIONS #
######################

cols_to_keep <- c("REACTANT_Name",
                  "PRODUCT_Name",
                  'Reaction.type',
                  'CFTR_interactions')
good_colnames <- c('Source',
                   'Target',
                   'Effect',
                   'CFTR_interactions')

WT_CFTR_interactions.PPI.connected.min <- WT_CFTR_interactions.PPI.connected[,cols_to_keep]
colnames(WT_CFTR_interactions.PPI.connected.min) <- good_colnames
### Remove duplicates
WT_CFTR_interactions.PPI.connected.min <- WT_CFTR_interactions.PPI.connected.min[!duplicated(WT_CFTR_interactions.PPI.connected.min),]


F508del_CFTR_interactions.PPI.connected.min <- F508del_CFTR_interactions.PPI.connected[,cols_to_keep]
colnames(F508del_CFTR_interactions.PPI.connected.min) <- good_colnames

### Replace complex names by "CFTR"
F508del_CFTR_interactions.PPI.connected.min[which(F508del_CFTR_interactions.PPI.connected.min$Source %in% F508del_CFTR), "Source"]<- "CFTR"
F508del_CFTR_interactions.PPI.connected.min[which(F508del_CFTR_interactions.PPI.connected.min$Target %in% F508del_CFTR), "Target"]<- "CFTR"
### Remove duplicates
F508del_CFTR_interactions.PPI.connected.min <- F508del_CFTR_interactions.PPI.connected.min[!duplicated(F508del_CFTR_interactions.PPI.connected.min),]

### To track the maps to which the interactions have been retrieved
WT_CFTR_interactions.PPI.connected.min$map <- "WT"
F508del_CFTR_interactions.PPI.connected.min$map <- "F508del"

### merge
all_CFTR_interactions.PPI.connected.min <- merge(WT_CFTR_interactions.PPI.connected.min,
                                                 F508del_CFTR_interactions.PPI.connected.min,
                                                 by = c("Source","Target","Effect", "CFTR_interactions"),   
                                                 all = TRUE)

### status 
all_CFTR_interactions.PPI.connected.min$status <- apply(all_CFTR_interactions.PPI.connected.min[,c("map.x", "map.y")],
                                                        MARGIN = 1,
                                                        FUN = function(x){
                                                          if (all(!is.na(x))){
                                                            return('unchanged')
                                                          } else if (is.na(x[1])){
                                                            return('gained')
                                                          } else {
                                                            return('lost')
                                                          }
                                                        })
all_CFTR_interactions.PPI.connected.min[,c("map.x", "map.y")] <- NULL
all_CFTR_interactions.PPI.connected.min$Effect <- as.character(all_CFTR_interactions.PPI.connected.min$Effect)

save(all_CFTR_interactions.PPI.connected.min,
     file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_interactions_df_2023_07_07.RData")

#########
# NODES #
#########

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

### Getting UniProt ID 
all_CFTR_interactions.PPI.connected.nodes.UniProtID <- getBM(filters= "hgnc_symbol", 
                                                             attributes= c("ensembl_gene_id","uniprotswissprot", "hgnc_symbol"),
                                                             values=unique(c(all_CFTR_interactions.PPI.connected.min$Source,
                                                                             all_CFTR_interactions.PPI.connected.min$Target)),
                                                             mart= mart)
all_CFTR_interactions.PPI.connected.nodes.UniProtID.unique <- all_CFTR_interactions.PPI.connected.nodes.UniProtID %>%
  filter(uniprotswissprot != "")
all_CFTR_interactions.PPI.connected.nodes.df <- unique(all_CFTR_interactions.PPI.connected.nodes.UniProtID.unique[,c("hgnc_symbol",
                                                                                                        "uniprotswissprot")])
colnames(all_CFTR_interactions.PPI.connected.nodes.df) <- c("HGNC", "UniProtID")

### CFTR interactors
all_CFTR_interactions.PPI <- all_CFTR_interactions.PPI.connected.min[which(all_CFTR_interactions.PPI.connected.min$CFTR_interactions),]
all_CFTR_interactors.PPI <- setdiff(unique(c(all_CFTR_interactions.PPI$Source,
                                             all_CFTR_interactions.PPI$Target)),
                                    c("CFTR"))
all_CFTR_interactions.PPI.connected.nodes.df$CFTR_interactor <-
  all_CFTR_interactions.PPI.connected.nodes.df$HGNC %in% all_CFTR_interactors.PPI

save(all_CFTR_interactions.PPI.connected.nodes.df,
     file = "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/diff_pathways_networks/CFTR_interactors_nodes_df_2023_07_07.RData")
