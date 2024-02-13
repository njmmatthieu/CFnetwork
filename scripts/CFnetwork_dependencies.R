if (!requireNamespace("data.table")){
  install.packages("BiocManager")
}

if (!requireNamespace("dplyr")){
  install.packages("dplyr")
}

if (!requireNamespace("igraph")){
  install.packages("igraph")
}

if (!requireNamespace("ggplot2")){
  install.packages("ggplot2")
}

if (!requireNamespace("pheatmap")){
  install.packages("pheatmap")
}

if (!requireNamespace("RColorBrewer")){
  install.packages("RColorBrewer")
}

if (!requireNamespace("tibble")){
  install.packages("tibble")
}

if (!requireNamespace("tidyr")){
  install.packages("tidyr")
}

if (!requireNamespace("tidyverse")){
  install.packages("tidyverse")
}

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install("biomaRt")
BiocManager::install("dorothea")
BiocManager::install("fgsea")
BiocManager::install("OmnipathR")
