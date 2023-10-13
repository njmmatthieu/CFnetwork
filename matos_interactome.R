library("biomaRt")
library(readxl)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

Matos_sup_mat_filename <- "/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/comparison_RegActinCyto/Matos_2019_sup_mat/mmc1/144713_2_supp_356136_ptymfq.xlsx"

excel_sheets("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/comparison_RegActinCyto/Matos_2019_sup_mat/mmc1/144713_2_supp_356136_ptymfq.xlsx")

wt_CFTR_interactome.df <- read_excel(Matos_sup_mat_filename,
                                     sheet = "Table S6",skip = 2)
wt_CFTR_interactome.HGNC.Df <- getBM(filters= "uniprotswissprot", 
                   attributes= c("ensembl_gene_id","uniprotswissprot", "hgnc_symbol"),
                   values=wt_CFTR_interactome.df$Uniprot_ID,
                   mart= mart)
wt_CFTR_interactome.Symbol <- unique(wt_CFTR_interactome.HGNC.Df$hgnc_symbol) 

F508del_interactome.df <- read_excel(Matos_sup_mat_filename,
                                     sheet = "Table S7",skip = 2)
F508del_interactome.HGNC.Df <- getBM(filters= "uniprotswissprot", 
                                     attributes= c("ensembl_gene_id","uniprotswissprot", "hgnc_symbol"),
                                     values=F508del_interactome.df$Uniprot_ID,
                                     mart= mart)
F508del_interactome.Symbol <- unique(F508del_interactome.HGNC.Df$hgnc_symbol) 
