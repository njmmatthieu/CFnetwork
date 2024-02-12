Verhaeghe_deg_df <- read.table("differential_expression_data/Verhaeghe_DEG_limma_unique_HGNC.txt",
                               sep = "\t", 
                               header = T)

Voisin_deg_df <- read.table("differential_expression_data/Voisin_DEG_limma_unique_HGNC.txt",
                            sep = "\t", 
                            header = T)

Ogilvie_nasal_deg_df <- read.table("differential_expression_data/Ogilvie_nasal_DEG_limma_unique_HGNC.csv",
                                   sep = "\t", 
                                   header = T)

Clarke_deg_df <- read.table("differential_expression_data/Clarke_DEG_limma_unique_HGNC.csv",
                            sep = "\t", 
                            header = T)

Ogilvie_bronchial_deg_df <- read.table("differential_expression_data/Ogilvie_bronchial_DEG_limma_unique_HGNC.csv",
                                       sep = "\t", 
                                       header = T)

Ogilvie_bronchial_deg_df$diffexpressed <- "NO"
Ogilvie_bronchial_deg_df$diffexpressed[Ogilvie_bronchial_deg_df$logFC > 1.0 & Ogilvie_bronchial_deg_df$adj.P.Val < 0.05] <- "UP"
Ogilvie_bronchial_deg_df$diffexpressed[Ogilvie_bronchial_deg_df$logFC < -1.0 & Ogilvie_bronchial_deg_df$adj.P.Val < 0.05] <- "DOWN"

Balloy_deg_df <- read.table("differential_expression_data/Balloy_without_CTRL2_DEG_limma_unique_HGNC.txt",
                            sep = "\t", 
                            header = T)

Zoso_deg_df <- read.table("differential_expression_data/Zoso_DEG_limma_unique_HGNC.txt", 
                          sep = "\t", 
                          header = T)
Zoso_deg_df$diffexpressed <- "NO"
Zoso_deg_df$diffexpressed[Zoso_deg_df$logFC > 1.0 & Zoso_deg_df$adj.P.Val < 0.05] <- "UP"
Zoso_deg_df$diffexpressed[Zoso_deg_df$logFC < -1.0 & Zoso_deg_df$adj.P.Val < 0.05] <- "DOWN"

Ling_deg_df <- read.table("differential_expression_data/Ling_normalised_differential_expression_analysis_filtered_limma.txt",
                          sep = "\t", 
                          header = T)
Ling_deg_df$diffexpressed <- "NO"
Ling_deg_df$diffexpressed[Ling_deg_df$logFC > 1.0 & Ling_deg_df$adj.P.Val < 0.05] <- "UP"
Ling_deg_df$diffexpressed[Ling_deg_df$logFC < -1.0 & Ling_deg_df$adj.P.Val < 0.05] <- "DOWN"

Saint_Criq_SC_deg_df <- read.table("differential_expression_data/Saint_Criq_SC_DEG_limma_unique_HGNC.txt",
                                   sep = "\t", 
                                   header = T)
Saint_Criq_UNC_deg_df <- read.table("differential_expression_data/Saint_Criq_UNC_without_NCF_2_DEG_limma_unique_HGNC.txt",
                                    sep = "\t", 
                                    header = T)

deg_results_list <- list(Verhaeghe_deg_df,
                         Voisin_deg_df,
                         Clarke_deg_df,
                         Ogilvie_nasal_deg_df,
                         Ogilvie_bronchial_deg_df,
                         Balloy_deg_df,
                         Zoso_deg_df,
                         Ling_deg_df,
                         Saint_Criq_UNC_deg_df,
                         Saint_Criq_SC_deg_df)

studies_names <- c("Verhaeghe",
                   "Voisin",
                   "Clarke",
                   "Ogilvie nasal",
                   "Ogilvie bronchial",
                   "Balloy",
                   "Zoso",
                   "Ling",
                   "Saint-Criq UNC",
                   "Saint-Criq SC")

names(deg_results_list) <- studies_names

# rm(Verhaeghe_deg_df,
#      Voisin_deg_df,
#      Clarke_deg_df,
#      Ogilvie_nasal_deg_df,
#      Ogilvie_bronchial_deg_df,
#      Balloy_deg_df,
#      Zoso_deg_df,
#      Ling_deg_df,
#      Saint_Criq_SC_deg_df,
#      Saint_Criq_UNC_deg_df)

diff_expressed_colmun <- function(deg_df, logFC_thr=1, padj_thr=0.05){
  
  deg_output_df <- deg_df
  deg_output_df$diffexpressed <- "NO"
  deg_output_df$diffexpressed[deg_output_df$logFC > logFC_thr & deg_output_df$adj.P.Val < padj_thr] <- "UP"
  deg_output_df$diffexpressed[deg_output_df$logFC < -logFC_thr & deg_output_df$adj.P.Val < padj_thr] <- "DOWN"
  
  return(deg_output_df)
}

