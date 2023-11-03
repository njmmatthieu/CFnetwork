# Construction of the CF signaling network

This repository contains all the scripts and data to reproduce the construction and the analysis of the CF signalling network. The overall network can be accessed as a Cytoscape session, in the [sysbio-curie/CFnetwork_cystoscape](https://github.com/sysbio-curie/CFnetwork_cytoscape.git) github repository for further analysis. 

Scripts were designed by Matthieu Najm (matthieu dot najm at minesparis dot psl dot eu). See a complete list of authors in the corresponding paper.
Feel free to explore, run, and adapt these scripts to suit your analysis needs. If you have any questions or need assistance, please don't hesitate to reach out.

<img src="./images/graphical_abstract_2023_10_09.png" width="600">



Explore the following sections to understand the content and purpose of each file:

## 0. Databases Preprocessing

### 0.1 KEGG Pathways Preprocessing
#### Script
- **scripts/kegg_pathways_scripts/kegg_pathways_preprocess.R** - Retrieve KEGG signaling pathways.

#### Dependencies
- scripts/kegg_pathways_scripts/kegg_pathways_utils.R** - Download KEGG pathway data.

#### Outputs
- kegg_pathways/Kegg_from_omnipathR_gsea_2022_09_07.gmt
- kegg_pathways/kegg_pathwaysKegg_pathways_from_omnipathR.txt
- kegg_pathways/kegg_pathways_from_omnipath_list.Rdata
- kegg_pathways/Kegg_pathways_from_omnipath_nodes_carac.Rdata

#### Visualisation
- kegg_pathways/example_kegg_pathway.R

### 0.2 CyFi-MAP Preprocessing
#### Script
- **scripts/CyFi-MAP_scripts/extract_CyFi_MAP_CFTR_interactors.R** - Extract CFTR interactors from the CyFi database

#### Inputs
- CyFi-MAP/CFTR_cp17-elementExport_2022_10_26_MN_2.txt : elements of the wt-CFTR interactome 
- CyFi-MAP/CFTR_cp17-networkExport_2022_10_26.txt: network interactions of the wt-CFTR interactome
- CyFi-MAP/F508del_cp21-elementExport_MN.txt : elements of the F508del-CFTR interactome 
- CyFi-MAP/F508del_cp21-networkExport.txt : network interactions of the F508del-CFTR interactome

#### Dependencies
- scripts/CyFi-MAP_scripts/CyFi_MAP_helper.R : Transform complexes to proteins.

#### Outputs
- CFTR_interactors/CFTR_interactors_nodes_df_2023_07_07.RData : CFTR interactors extracted from the CyFi-MAP database
- CFTR_interactors/CFTR_interactors_interactions_df_2023_07_07.RData : CFTR interactions extracted from the CyFi-MAP database

## 1. Transcriptomics analyses 

### 1.1 Differential Expression Analysis at the pathway Level with fGSEA

#### Scripts
- scripts/fgsea/Verhaeghe_fgsea_2022_01_24.Rmd
- scripts/fgsea/Ogilvie_nasal_fgsea_2022_01_24.Rmd
- scripts/fgsea/Ogilvie_bronchial_fgsea_2022_01_14.Rmd
- scripts/fgsea/Voisin_fgsea_2022_01_19.Rmd
- scripts/fgsea/Clarke_fgsea_2022_01_14.Rmd
- scripts/fgsea/Balloy_fgsea_2022_02_25.Rmd
- scripts/fgsea/Zoso_fgsea_2022_01_14.Rmd
- scripts/fgsea/Ling_fgsea_2022_01_27.Rmd
- scripts/fgsea/Saint_Criq_fgsea_2022_01_27.Rmd

#### Inputs
##### Gene level t-statistics:
- differential_expression_data/Verhaeghe_DEG_limma_unique_HGNC_2022_02_07.txt
- differential_expression_data/Ogilvie_nasal_DEG_limma_unique_HGNC_2022_02_09.csv
- differential_expression_data/Ogilvie_bronchial_DEG_limma_unique_HGNC_2022_02_09.csv
- differential_expression_data/Voisin_DEG_limma_unique_HGNC_2021_02_08.txt
- differential_expression_data/Clarke_DEG_limma_unique_HGNC_2022_01_24.csv
- differential_expression_data/Balloy_without_CTRL2_DEG_limma_unique_HGNC_2022_02_25.txt
- differential_expression_data/Zoso_DEG_limma_unique_HGNC_2021_02_09.txt
- differential_expression_data/Saint_Criq_UNC_without_NCF_2_DEG_limma_unique_HGNC_2021_01_27.txt
- differential_expression_data/Saint_Criq_SC_DEG_limma_unique_HGNC_2021_01_27.txt

##### GMT:
- kegg_pathways/Kegg_from_omnipathR_gsea_2022_09_07.gmt

#### Outputs
- fgsea_output/Verhaeghe_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Ogilvie_nasal_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Ogilvie_bronchial_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Voisin_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Clarke_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Balloy_without_CTRL2_unfiltered_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Zoso_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Ling_unfiltered_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Saint_Criq_UNC_without_NCF_2_fgsea_report_kegg_pathways_2022_09_07.RData
- fgsea_output/Saint_Criq_SC_fgsea_report_kegg_pathways_2022_09_07.RData

### 1.2 Meta-analysis at the pathway level

#### Script
**scripts/fgsea_comparison/fGSEA_comparison.Rmd**: Perform pathway-level differential expression analysis (starting from line 282).

##### Dependencies: 
    - scripts/fgsea_comparison/Fgsea_output_preprocess.R
    - scripts/fgsea_comparison/Fgsea_comparison_utils.R

### Output
- fgsea_nes_diff_pathways_2022_09_15.RData

## 2. From pathways to CF Network

### Main Script
**kegg_diff_pathways.Rmd** : Correct and merge dysregulated pathways into a single network.

### Input: 
- kegg_pathways/kegg_diff_pathways_corrections_w_EZR_2023_07_07.txt: list of interactions and nodes to correct in KEGG pathways

#### Dependencies:
##### Correction of KEGG pathways:
- pathways_to_network/kegg_pathways_manual_curation.R: Correct KEGG signaling pathways.
- pathways_to_network/network_utils.R

##### Add CFTR interactors (Optional):
- pathways_to_network/CFTR_interactors_helper.R

##### Outputs:
- kegg_diff_pathway_network/kegg_diff_pathways_interactions_with_CFTR_interactors_df_2023_07_10.RData : CF network interactions as data.frame
- kegg_diff_pathway_network/kegg_diff_pathways_nodes_with_CFTR_interactors_df_2023_07_10.RData : CF network nodes as data.frame

## 3. CF network pruning

### Main Script:
**scripts/kegg_diff_pathways_network_scripts/kegg_diff_pathways_layout.R**

#### Input: 
- kegg_diff_pathway_network/kegg_diff_pathways_interactions_with_CFTR_interactors_df_2023_07_10.RData : CF network interactions as data.frame
- kegg_diff_pathway_network/kegg_diff_pathways_nodes_with_CFTR_interactors_df_2023_07_10.RData : CF network nodes as data.frame

#### Dependencies:
- scripts/kegg_diff_pathways_network_scripts/simplify_netowrk_helper.R 
- scripts/kegg_diff_pathways_network_scripts/network_visualisation_helper.R

#### Outputs:
- kegg_diff_pathway_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df_2023_07_10.txt : final pruned CF network interactions as txt file.
- kegg_diff_pathway_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df_2022_07_10.txt : final pruned CF network nodes as txt file.

## 4. CF network analysis

### Main Script:
**kegg_diff_pathways_network_analysis_final.R** - Network analysis 

#### Input: 
- kegg_diff_pathway_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df_2023_07_10.txt : final pruned CF network interactions as txt file.
- kegg_diff_pathway_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df_2022_07_10.txt : final pruned CF network nodes as txt file.




## A. Supplementary Analysis

## A.1. - Exploratory analysis

## A.2. -  Differential Expression Analysis at the gene Level

### Main Script
**scripts/DEG_comparison/DEG_comparison.Rmd** : Compare differentially expressed genes.

#### Dependencies: 
- scripts/DEG_comparison/deg_utils.R







