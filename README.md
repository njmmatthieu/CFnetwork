# Construction of the CF signaling network

This repository contains all the scripts and data to reproduce the construction and the analysis of the CF signalling network. The overall network can be accessed as a Cytoscape session, in the [sysbio-curie/CFnetwork_cystoscape](https://github.com/sysbio-curie/CFnetwork_cytoscape.git) github repository for further analysis. 

Please note that the repository is currently undergoing some cleanup and organization. It may not be in its final, clean state. I'm currently working on making this repository more organized and user-friendly. If you have any questions or suggestions, please don't hesitate to reach out (matthieu dot najm at minesparis dot psl dot eu).


Feel free to explore, run, and adapt these scripts to suit your analysis needs.
Explore the following sections to understand the content and purpose of each file:

## 1. Transcriptomics analyses 

### 1.1 Differential Expression Analysis at the pathway level with fGSEA

Perform pathway-level differential expression analysis with fgsea algorithm for each study.

#### Conda environment
fgsea

#### Scripts
    - scripts/fgsea/Verhaeghe_fgsea.Rmd
    - scripts/fgsea/Ogilvie_nasal_fgsea.Rmd
    - scripts/fgsea/Ogilvie_bronchial.Rmd
    - scripts/fgsea/Voisin_fgsea.Rmd
    - scripts/fgsea/Clarke_fgsea.Rmd
    - scripts/fgsea/Balloy_fgsea.Rmd
    - scripts/fgsea/Zoso_fgsea.Rmd
    - scripts/fgsea/Ling_fgsea.Rmd
    - scripts/fgsea/Saint_Criq_fgsea.Rmd

#### Inputs
##### Gene level t-statistics:
    - data/differential_expression_data/Verhaeghe_DEG_limma_unique_HGNC.txt
    - data/differential_expression_data/Ogilvie_nasal_DEG_limma_unique_HGNC.csv
    - data/differential_expression_data/Ogilvie_bronchial_DEG_limma_unique_HGNC.csv
    - data/differential_expression_data/Voisin_DEG_limma_unique_HGNC.txt
    - data/differential_expression_data/Clarke_DEG_limma_unique_HGNC.csv
    - data/differential_expression_data/Balloy_without_CTRL2_DEG_limma_unique_HGNC.txt
    - data/differential_expression_data/Zoso_DEG_limma_unique_HGNC.txt
    - data/differential_expression_data/Saint_Criq_UNC_without_NCF_2_DEG_limma_unique_HGNC.txt
    - data/differential_expression_data/Saint_Criq_SC_DEG_limma_unique_HGNC.txt

##### GMT:
    - data/kegg_pathways/kegg_pathways_from_omnipathR.gmt

#### Outputs
    - data/fgsea_output/Verhaeghe_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Ogilvie_nasal_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Ogilvie_bronchial_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Voisin_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Clarke_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Balloy_without_CTRL2_unfiltered_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Zoso_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Ling_unfiltered_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Saint_Criq_UNC_without_NCF_2_fgsea_report_kegg_pathways.RData
    - data/fgsea_output/Saint_Criq_SC_fgsea_report_kegg_pathways.RData

### 1.2 Meta-analysis at the pathway level

Compare the results of the pathway level analysis of all the studies to extract a list of common dysregulated pathways in CF.

#### Conda environment
deseq

#### Script
    - scripts/fgsea_comparison/fGSEA_comparison.Rmd

#### Dependencies 
    - scripts/fgsea_comparison/Fgsea_output_preprocess.R
    - scripts/fgsea_comparison/Fgsea_comparison_utils.R

#### Output
    - data/kegg_diff_pathways/fgsea_nes_diff_pathways.RData

## 2. From pathways to CF Network

Correct and merge the common dysregulated pathways into one single network. Add CFTR and CFTR interactors to the network.

#### Conda environment
deseq

#### Script
    - scripts/pathways_to_network/kegg_diff_pathways_to_network.Rmd

#### Input 
List of interactions and nodes to correct in KEGG pathways
    - data/kegg_pathways/kegg_diff_pathways_corrections_w_EZR_2023_07_07.txt

#### Dependencies
##### Correction of KEGG pathways
    - scripts/kegg_pathways/kegg_pathways_manual_curation.R
    - scripts/pathways_to_network/network_utils.R

##### Add CFTR interactors (Optional)
    - scripts/pathways_to_network/CFTR_interactors_helper.R

##### Outputs
    - data/kegg_diff_pathways_network/kegg_diff_pathways_interactions_with_CFTR_interactors_df.RData : CF network interactions as data.frame
    - data/kegg_diff_pathways_network/kegg_diff_pathways_nodes_with_CFTR_interactors_df.RData : CF network nodes as data.frame

## 3. CF network pruning

Network pruning and cleaning (corresponding to the part 3.6 of the Methods section). 

#### Script
    - scripts/network_analysis/kegg_diff_pathways_network_layout.R

#### Input
    - data/kegg_diff_pathways_network/kegg_diff_pathways_interactions_with_CFTR_interactors_df.RData : CF network interactions as data.frame
    - data/kegg_diff_pathways_network/kegg_diff_pathways_nodes_with_CFTR_interactors_df.RData : CF network nodes as data.frame

#### Dependencies
    - scripts/network_analysis/simplify_network_helper.R 
    - scripts/network_analysis/network_visualisation_helper.R

#### Outputs
    - data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df.txt : final pruned CF network interactions as txt file.
    - data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df.txt : final pruned CF network nodes as txt file.

## 4. CF network analysis

Network topological analysis (corresponding to the part 2.5 of the Results). 

#### Script
    - scripts/network_analysis/kegg_diff_pathways_network_analysis_final.R

#### Input 
    - data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_interactions_df.txt : final pruned CF network interactions as txt file.
    - data/kegg_diff_pathways_network/diff_kegg_pathways_with_CFTR_interactors_PPI_direct_tagged_nodes_df.txt : final pruned CF network nodes as txt file.

\
\
\

## A. Supplementary Analysis

### A.1 KEGG Pathways Preprocessing

#### Conda environment
doro-r4

#### Script - Retrieve KEGG signaling pathways.
    - scripts/kegg_pathways/kegg_pathways_preprocess.R

#### Outputs
    - kegg_pathways/kegg_pathways_from_omnipathR.gmt
    - kegg_pathways/kegg_pathwaysKegg_pathways_from_omnipathR.txt
    - kegg_pathways/kegg_pathways_from_omnipath_list.Rdata
    - kegg_pathways/Kegg_pathways_from_omnipath_nodes_carac.Rdata

#### Visualisation
    - kegg_pathways/example_kegg_pathway.R

### A.2 CyFi-MAP Preprocessing

#### Conda environment
deseq

#### Script - Extract CFTR interactors from the CyFi database
    - scripts/CyFi-MAP_scripts/extract_CyFi_MAP_CFTR_interactors.R 

#### Inputs
    - CyFi-MAP/CFTR_cp17-elementExport_2022_10_26_MN_2.txt : elements of the wt-CFTR interactome. 
    - CyFi-MAP/CFTR_cp17-networkExport_2022_10_26.txt: network interactions of the wt-CFTR interactome.
    - CyFi-MAP/F508del_cp21-elementExport_MN.txt : elements of the F508del-CFTR interactome. 
    - CyFi-MAP/F508del_cp21-networkExport.txt : network interactions of the F508del-CFTR interactome.

#### Dependencies - Transform complexes to proteins.
    - scripts/CyFi-MAP_scripts/CyFi_MAP_helper.R

#### Outputs
    - CFTR_interactors/CFTR_interactors_nodes_df.RData : CFTR interactors extracted from the CyFi-MAP database.
    - CFTR_interactors/CFTR_interactors_interactions_df.RData : CFTR interactions extracted from the CyFi-MAP database.

### A.3. -  Differential Expression Analysis at the gene Level

#### Script - Compare differentially expressed genes.
    - scripts/DEG_comparison/DEG_comparison.Rmd

#### Dependencies 
    - scripts/DEG_comparison/deg_utils.R







