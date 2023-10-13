library(cluster)
library(pheatmap)
library(tidyverse)

# GSEA

load("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/diff_kegg_pathways_nodes_df_2022_09_08.RData")

kegg_diff_pathways_nodes_wide.only_pathway <- kegg_diff_pathways_nodes_wide[,!(colnames(kegg_diff_pathways_nodes_wide) %in% c("sum",
                                                                                "pathway"))]
rownames(kegg_diff_pathways_nodes_wide.only_pathway) <- kegg_diff_pathways_nodes_wide.only_pathway$Symbol
kegg_diff_pathways_nodes_wide.only_pathway$Symbol <- NULL
GSEA <- kegg_diff_pathways_nodes_wide.only_pathway



# Leading Edges

leadingEdges_nodes <- read.table("/Users/matthieu/ownCloud/Thèse/Systems Biology/pathways_to_network/networks/pathway_networks/kegg_pathways/diff_kegg_pathways_leadingEdges_nodes_2022_09_15.txt",
                                 sep = "\t",
                                 header = T)

leadingEdges_nodes.only_pathway <- leadingEdges_nodes[, !(colnames(leadingEdges_nodes) %in% c("sum",
                                                                                             "pathway",
                                                                                             "tf_targets",
                                                                                             "Verhaeghe",
                                                                                             "Voisin",
                                                                                             "Clarke",
                                                                                             "Ogilvie.nasal",
                                                                                             "Ogilvie.bronchial",
                                                                                             "Balloy",
                                                                                             "Zoso",
                                                                                             "Ling",
                                                                                             "Saint.Criq.UNC",
                                                                                             "Saint.Criq.SC"))]
rownames(leadingEdges_nodes.only_pathway) <- leadingEdges_nodes.only_pathway$Symbol
leadingEdges_nodes.only_pathway$Symbol <- NULL
GSEA <- leadingEdges_nodes.only_pathway

# clustering and visualisatiopn

pathways_dist_mat <- dist(t(GSEA), 
                          method = 'binary')
pathways_clust_avg <- hclust(pathways_dist_mat, 
                             method = 'average')

genes_dist_mat <- dist(GSEA, 
                          method = 'binary')
genes_clust_avg <- hclust(genes_dist_mat, 
                             method = 'average')

leadingEdges_nodes.ordered_pathway <- GSEA[genes_clust_avg$labels[genes_clust_avg$order],
                                           pathways_clust_avg$labels[pathways_clust_avg$order]]

pheatmap(leadingEdges_nodes.ordered_pathway,
         cluster_rows = F,
         cluster_cols = F)
distfunc <- function(x) daisy(x, metric = "euclidean")
d <- distfunc(leadingEdges_nodes.only_pathway)

clust <- pheatmap(as.matrix(d),
                  show_rownames = T)

nb_common_genes_all_pathways_df <- data.frame()
for (pathway_A in colnames(GSEA)){
  # print(pathway_A)
  for (pathway_B in colnames(GSEA)){
    # print(pathway_B)
    if (pathway_A==pathway_B){
      nb_common_genes <- NA
    } else {
      nb_common_genes <- sum(as.integer(sapply(1:nrow(GSEA), function(i_row) {return(all(GSEA[i_row,c(pathway_A, pathway_B)]==1))})))
    }
    nb_common_genes_df <- data.frame(pathway_A = pathway_A,
                                     pathway_B = pathway_B,
                                     nb_common_genes = nb_common_genes)
    # print(nb_common_genes_df)
    nb_common_genes_all_pathways_df <- rbind(nb_common_genes_all_pathways_df,
                                             nb_common_genes_df)
    
}
}

nb_common_genes_all_pathways_wide <- nb_common_genes_all_pathways_df %>% 
  spread(pathway_B, nb_common_genes)
rownames(nb_common_genes_all_pathways_wide) <- nb_common_genes_all_pathways_wide$pathway_A
nb_common_genes_all_pathways_wide$pathway_A <- NULL


nb_common_genes_all_pathways_wide$mean_common_genes <- sapply(1:nrow(nb_common_genes_all_pathways_wide), function(i_row){
  return(mean(as.numeric(nb_common_genes_all_pathways_wide[i_row,]), na.rm=T))
})

nb_common_genes_all_pathways_wide$nb_non_connected_pathways <- sapply(1:nrow(nb_common_genes_all_pathways_wide), function(i_row){
  return(sum(as.integer(nb_common_genes_all_pathways_wide[i_row,]==0), na.rm=T))
})

nb_common_genes_all_pathways_wide$nb_less_connected_pathways <- sapply(1:nrow(nb_common_genes_all_pathways_wide), function(i_row){
  return(sum(as.integer(nb_common_genes_all_pathways_wide[i_row,]<=1), na.rm=T))
})


pathways_clust_avg <- hclust(nb_common_genes_all_pathways_wide, 
                             method = 'average')

pheatmap(nb_common_genes_all_pathways_wide,
         cluster_rows = T)

pathway_A <- colnames(GSEA)[1]
pathway_B <- colnames(GSEA)[2]
