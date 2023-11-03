
# A - Micro-Arrays Cell Lines 

## 1 - Verhaeghe

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Verhaeghe/GSEA/Verhaeghe_fgsea_report_kegg_pathways_2022_09_07.RData")
Verhaeghe_fgsea <- fgseaRes 

## 2 - Voisin

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Voisin/GSEA/Voisin_fgsea_report_kegg_pathways_2022_09_07.RData")
Voisin_fgsea <- fgseaRes 

# B - Micro-array nasal - Primary cultures

## 3 - Ogilvie nasal

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Ogilvie/GSEA/Ogilvie_nasal_fgsea_report_kegg_pathways_2022_09_07.RData")
Ogilvie_nasal_fgsea <- fgseaRes 

## 4 - Clarke

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Clarke/GSEA/Clarke_fgsea_report_kegg_pathways_2022_09_07.RData")
Clarke_fgsea <- fgseaRes 

# C - Micro-array bronchial - Primary cultures

## 5 - Ogilvie bronchial

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Ogilvie/GSEA/Ogilvie_bronchial_fgsea_report_kegg_pathways_2022_09_07.RData")
Ogilvie_bronchial_fgsea <- fgseaRes 

# D - RNA seq - bronchial

## 6 - Balloy

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Balloy/GSEA/Balloy_without_CTRL2_unfiltered_fgsea_report_kegg_pathways_2022_09_07.RData")
Balloy_fgsea <- fgseaRes 

## 7 - Zoso

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Zoso/GSEA/Zoso_fgsea_report_kegg_pathways_2022_09_07.RData")
Zoso_fgsea <- Zoso_limma_fgseaRes 


## 8 - Ling

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Ling/GSEA/Ling_unfiltered_fgsea_report_kegg_pathways_2022_09_07.RData")
Ling_fgsea <- fgseaRes 

## 9 - Saint-Criq (SC)

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Saint-Criq/GSEA/Saint_Criq_SC_fgsea_report_kegg_pathways_2022_09_07.RData")
Saint_Criq_SC_fgsea <- fgseaRes 


## 10 - Saint-Criq (UNC)

load("~/ownCloud/Thèse/Systems Biology/Transcriptomic studies/Transcriptomic Analyses/Analyse Saint-Criq/GSEA/Saint_Criq_UNC_without_NCF_2_fgsea_report_kegg_pathways_2022_09_07.RData")
Saint_Criq_UNC_fgsea <- fgseaRes 

# fgsea_results_list <- list(Verhaeghe_fgsea[,-which(colnames(Verhaeghe_fgsea) %in% c("leadingEdge"))],
#                            Voisin_fgsea[,-which(colnames(Voisin_fgsea) %in% c("leadingEdge"))],
#                            Ogilvie_nasal_fgsea[,-which(colnames(Ogilvie_nasal_fgsea) %in% c("leadingEdge"))],
#                            Clarke_fgsea[,-which(colnames(Clarke_fgsea) %in% c("leadingEdge"))],
#                            Ogilvie_bronchial_fgsea[,-which(colnames(Ogilvie_bronchial_fgsea) %in% c("leadingEdge"))],
#                            Balloy_fgsea[,-which(colnames(Balloy_fgsea) %in% c("leadingEdge"))],
#                            Zoso_fgsea[,-which(colnames(Zoso_fgsea) %in% c("leadingEdge"))],
#                            Ling_fgsea[,-which(colnames(Ling_fgsea) %in% c("leadingEdge"))],
#                            Saint_Criq_UNC_fgsea[,-which(colnames(Saint_Criq_UNC_fgsea) %in% c("leadingEdge"))],
#                            Saint_Criq_SC_fgsea[,-which(colnames(Saint_Criq_SC_fgsea) %in% c("leadingEdge"))])




# fgsea_results_list <- lapply(fgsea_results_list, function(fgsea_results_list){
#   
# })

fgsea_diff_activated <- function(fgsea_output, padj_threshold=0.05) {
   
   fgsea_output$diff_activated<- "NO"
   fgsea_output$diff_activated[fgsea_output$NES > 0 & fgsea_output$padj < padj_threshold] <- "UP"
   fgsea_output$diff_activated[fgsea_output$NES < 0 & fgsea_output$padj < padj_threshold] <- "DOWN"
   
   return(fgsea_output)
   
}

fgsea_raw_results_list <- list(Verhaeghe_fgsea,
                               Voisin_fgsea,
                               Ogilvie_nasal_fgsea,
                               Clarke_fgsea,
                               Ogilvie_bronchial_fgsea,
                               Balloy_fgsea,
                               Zoso_fgsea,
                               Ling_fgsea,
                               Saint_Criq_UNC_fgsea,
                               Saint_Criq_SC_fgsea)

fgsea_results_list <- lapply(fgsea_raw_results_list, function(df){
  return(df[,-8])
}
)

studies_names <- c("Verhaeghe",
                   "Voisin",
                   "Ogilvie nasal",
                   "Ogilvie bronchial",
                   "Clarke",
                   "Balloy",
                   "Zoso",
                   "Ling",
                   "Saint-Criq UNC",
                   "Saint-Criq SC")

names(fgsea_results_list) <- studies_names

rm(fgseaRes,
   Zoso_limma_fgseaRes,
   Verhaeghe_fgsea,
     Voisin_fgsea,
     Ogilvie_nasal_fgsea,
     Clarke_fgsea,
     Ogilvie_bronchial_fgsea,
     Balloy_fgsea,
     Zoso_fgsea,
     Ling_fgsea,
     Saint_Criq_UNC_fgsea,
     Saint_Criq_SC_fgsea,
   fgsea_raw_results_list)
