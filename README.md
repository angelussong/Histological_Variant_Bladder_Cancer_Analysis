Some base R scripts using Seruat and monocle3 to analyze scRNA-seq data

Bladder_Seqwell: DGE matrix preprocessing, QC, clustering, DEG, and annotation

Epithelial_Seqwell: Epithelial/Tumor cell subsetting/sub-clustering, DEG, and annotation

InferCNV_Preparation: Make the matrix and phenotype file for running inferCNV

Correlation_Analysis.R: Correlation analysis of tumor cells with Cluster13

Cluster13_analysis: Handle the detailed analysis for cluster 13, including marker analyses and multiple ways of correlation analysis. 

monocle3_individual_tumor: pseudotime analysis for each variant tumor to assess the role of cluster13 in the pseudotime trajectory. 

Micropapillary_Nested_DEG: Derive common genes for micropapillary and nested variant tumors

Dividing_cell_analysis: Subset and characterize the dividing cells

GO_Analysis.R: API to run GSEA analysis for DEGs (with log2FC and p-value)

Specific_Luminal_Basal_CAM: Testing basal-luminal axis in the scRNA-seq data (not included in the paper but super interesting)
