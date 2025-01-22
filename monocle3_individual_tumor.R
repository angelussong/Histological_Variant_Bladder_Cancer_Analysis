library(monocle)
library(DESeq2)
library(edgeR)
library(DEFormats)
library(Seurat)
library(cowplot)
library(dplyr)
library(GEOquery)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ArrayTools)
library(dplyr)
library(monocle3)
# This script is designed to pipeline Seurat object output into Monocle3
# This script also contains all the necessary functionalities in Monocle3 as 
# This could also convert a 2D Seurat object and visualize/analyze it in 3D
# It will start from reading in a Seurat object with the count sparse matrix, UMAP coordinates, and cluster information
# This script is originally written for local machines but adaptations have also been included in annotations


### Require:: 'BiocManager', 'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats', 'limma', 'S4Vectors', 'SingleCellExperiment', 'SummarizedExperiment', 'reticulate', 'htmlwidgets'

# This is a required python package
#reticulate::py_install("louvain")

# This is installing the actual monocle3
#devtools::install_github('cole-trapnell-lab/monocle3')

library(jjb)
library(Seurat)
library(htmlwidgets)
library(monocle3)

### Installing the packages
setwd("~/Desktop/Bladder_2023/Pseudotime/")

dge<-readRDS("./dge_VAR08.rds")
mkdir("./VAR08")
setwd("./VAR08")



seurat<-dge
seurat$refinedID<-paste0("Cluster_",as.vector(seurat@active.ident))

#Extract data, phenotype data, and feature data from the SeuratObject
gene_annotation <- as.data.frame(rownames(seurat@reductions[["pca"]]@feature.loadings), row.names = rownames(seurat@reductions[["pca"]]@feature.loadings))
gene_annotation <- as.data.frame(rownames(seurat), row.names = rownames(seurat))

colnames(gene_annotation) <- "gene_short_name"

# part two, cell information

cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"
cell_metadata<-cbind(cell_metadata,as.vector(seurat$refinedID))
colnames(cell_metadata)<-c("barcode","refinedID")
# part three, counts sparse matrix

New_matrix <- seurat@assays[["RNA"]]@counts
#New_matrix <- seurat@assays[["RNA"]]@data
#New_matrix <- New_matrix[rownames(seurat@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix


### Construct the basic cds object

HSMM <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


# recreate.partition <- c(rep(1, length(HSMM@colData@rownames)))
# names(recreate.partition) <- HSMM@colData@rownames
# recreate.partition <- as.factor(recreate.partition)
# 
# HSMM@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info

# list_cluster <- seurat@meta.data[[sprintf("refinedID")]]
# #list_cluster <- seurat@meta.data[[sprintf("ClusterNames_%s_%sPC", 0.5, 20)]]
# names(list_cluster) <- seurat@assays[["RNA"]]@data@Dimnames[[2]]

# HSMM@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
# 
# 
# ### Could be a space-holder, but essentially fills out louvain parameters
# 
# HSMM@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"

#View data
#pData(HSMM)
#fData(HSMM)

HSMM <- estimate_size_factors(HSMM)

#print(dim(exprs(HSMM)))

## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM<-preprocess_cds(HSMM,method = "PCA")
#HSMM <- reduce_dimension(HSMM)

HSMM = reduce_dimension(HSMM, max_components = 2,reduction_method = "UMAP")
HSMM = cluster_cells(HSMM,reduction_method = "UMAP",resolution = 1e-3)
HSMM = learn_graph(HSMM,verbose = T,use_partition = F)


### Assign feature loading for downstream module analysis





# plot_cells(HSMM, 
#            color_cells_by = "partition")


p1<-plot_cells(HSMM,group_label_size = 3,cell_size = 1)
p2<-plot_cells(HSMM,color_cells_by = "refinedID", group_label_size = 3,cell_size = 1)
p3<-p1+p2
p3



ciliated_cds_pr_test_res <- graph_test(HSMM, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))


gene_module_df <- find_gene_modules(HSMM[pr_deg_ids,], resolution=c(0.0001,10^seq(-3,-1)))
colData(HSMM)$assigned_cell_type <- as.character(clusters(HSMM)[colnames(HSMM)])


cell_group_df <- tibble::tibble(cell=row.names(colData(HSMM)), 
                                cell_group=colData(HSMM)$refinedID)


agg_mat <- aggregate_gene_expression(HSMM, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
fig5<-pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")
pdf(paste0("Tumor_DE.pdf"), width = 15, height = 10)
fig5
dev.off()


# write.table(gene_module_df[order(gene_module_df$module,decreasing = F),],"DE_module.txt",col.names = T,row.names = T,quote = F,sep = "\t")
# module<-gene_module_df[order(gene_module_df$module,decreasing = F),]
# module$number<-as.vector(rownames(module))
# top10_module <- module %>% group_by(module) %>% top_n(n = -10,wt = number)
# write.table(top10_module,"DE_module_top10.txt",col.names = T,row.names = T,quote = F,sep = "\t")

mkdir("./Expression")
Genelist<-c("KRT20","TM4SF1","TOP2A","MKI67","KRT7","KRT24","MUC16","KRT8","KRT18","PTPRC","POU2F3","CD44")
Genelist<-c("TM4SF1","KRT24","WISP2","KRT7","KRT20","HLA-DRA","MUC16",
            "CHGA","KRT8","KRT19","KRT14")
for (i in 1:length(Genelist)){
  if(Genelist[i] %in% rownames(HSMM)){
  fig1<-plot_cells(HSMM,
                   genes=Genelist[i],
                   label_cell_groups=T,
                   show_trajectory_graph=T,cell_size = 1,min_expr = 0,scale_to_range = F)
  save_plot_pdf(fig1,paste0("./Expression/",Genelist[i],".pdf"),width = 6,height = 6)}
}

fig2<-plot_cells(HSMM,1,2,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,cell_size = 1)+facet_wrap(~refinedID)
pdf(paste0("Seurat_separate.pdf"), width = 15, height = 10)
fig2
dev.off()

HSMM = order_cells(HSMM)

# plot_cells(HSMM, 
#            color_cells_by = "pseudotime")

p1<-plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = F,label_groups_by_cluster = T,cell_size = 1, group_label_size = 3,show_trajectory_graph = F)
p2<-plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "pseudotime",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1, group_label_size = 3,show_trajectory_graph = F)
p3<-p1+p2
pdf("ID_with_Pseudotime.pdf", width = 10, height = 6)
p3
dev.off()


# p1<-plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "assigned_cell_type",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)
# p2<-plot_cells(HSMM,reduction_method = "UMAP",color_cells_by = "pseudotime",label_cell_groups = T,label_roots = T,label_groups_by_cluster = T,cell_size = 1)
# p3<-p1+p2
# pdf("Monocle_cluster_with_Pseudotime.pdf", width = 25, height = 10)
# p3
# dev.off()
# 
# 
# fig2<-plot_cells(HSMM,1,2,reduction_method = "UMAP",color_cells_by = "refinedID",label_cell_groups = T,label_roots = T,cell_size = 1, group_label_size = 5)+facet_wrap(~refinedID)
# pdf(paste0("Seurat_separate.pdf"), width = 30, height = 30)
# fig2
# dev.off()
# 
# fig2<-plot_cells(HSMM,1,2,reduction_method = "UMAP",color_cells_by = "pseudotime",label_cell_groups = T,label_roots = T,cell_size = 1, group_label_size = 5)+facet_wrap(~assigned_cell_type)
# pdf(paste0("Pseudotime_separate.pdf"), width = 30, height = 30)
# fig2
# dev.off()
# 
# genename<-c("MNP_a","MNP_b","MNP_c","MNP_d","Neutrophil")
# genes<-list()
# genes[[1]]<-c("S100A12","VCAN","CD163","THBS1","EREG")
# genes[[2]]<-c("TCF7L2","CDKN1C","CD79B","FCGR3A","RHOC")
# genes[[3]]<-c("FCER1A","RGS1","ADAM28","GSN","FCGR2B")
# genes[[4]]<-c("SEPP1","C1QC","C1QB","RNASE1","APOE")
# genes[[5]]<-c("FCGR3B","CMTM2","PTGS2","S100P","G0S2")
# 
# 
# 
# fig6<-plot_genes_in_pseudotime(lineage_HSMM,cell_size = 0.5,color_cells_by = "refinedID")
# pdf(paste0(genes,"_bygroup.pdf"), width = 6, height = 6)
# fig6
# # plot_genes_in_pseudotime(lineage_HSMM)
# dev.off()
# 
# 
# HSMM_metadata<-as.data.frame(colData(HSMM))
# Cells_in_question<-HSMM_metadata$barcode[(HSMM_metadata$assigned_cell_type=="2") & (HSMM_metadata$ID %in% c("Monocyte","MNPa"))]
# HSMM@colData$Size_Factor

Pseudotime<-pseudotime(HSMM, reduction_method = "UMAP")
seurat$Pseudotime<-Pseudotime

genelist<-VariableFeatures(seurat)
genelist<-genelist[!(genelist=="ID2")]
DF_cor<-NULL
DF_cor$Gene<-genelist
DF_cor$Corr<-0
DF_cor$pvalue<-1
DF_cor<-as.data.frame(DF_cor)

for (i in 1:length(genelist)){
  DF_tmp<-FetchData(seurat,vars = genelist[i])
  DF_tmp<-as.numeric(DF_tmp[,1])
  COR_test<-cor.test(DF_tmp,Pseudotime)
  DF_cor$Corr[i]<-COR_test$estimate
  DF_cor$pvalue[i]<-COR_test$p.value
}
DF_cor<-DF_cor[order(abs(DF_cor$Corr),decreasing = T),]
write.table(DF_cor,"./Correlation_Pseudotime.txt",sep="\t",row.names = F,col.names = T,quote = F)

mkdir("./Expression_bypseudotime")
Top50_cor<-DF_cor$Gene[1:50]
# for (i in 1:length(Top50_cor)){
# Gene=Top50_cor[i]
# Gene_expression<-FetchData(seurat,vars = Gene)
# Pseudotime_point<-Pseudotime
# Meta<-seurat$refinedID
# #cor.test(as.numeric(Expression_Gene1[,1]),as.numeric(Expression_Gene2[,1]))
# DF<-cbind(Gene_expression,Pseudotime_point,Meta)
# colnames(DF)<-c("Gene","Pseudotime","Group")
# DF$Gene<-as.numeric(DF$Gene)
# DF$Pseudotime<-as.numeric(DF$Pseudotime)
# ggplot(DF, aes(y = Pseudotime, x = Gene)) +
#   geom_point(aes(colour = factor(Group))) +
#   facet_wrap("Group") +
#   geom_smooth(method=lm) + 
#   xlab(Gene) + ylab("Pseudotime") +NoLegend()
# ggsave(paste0("./Expression_bypseudotime/Correlation_","Pseudotime","_",Gene,".pdf"),width = 6,height = 4)
# }


Genelist<-c("KRT7","CD44","KRT20","HOXB3","HOXB4","HOXB6","PRDM1","IL6R","SDC1")
for (i in 1:length(Genelist)){
  if(Genelist[i] %in% rownames(HSMM)){
    lineage_HSMM <- HSMM[rowData(HSMM)$gene_short_name %in% Genelist[i],]
    
    fig6<-plot_genes_in_pseudotime(lineage_HSMM,cell_size = 1,min_expr = 0,vertical_jitter = T,horizontal_jitter = T)
    
    save_plot_pdf(fig6,paste0("./Correlation_",Genelist[i],"_Linear.pdf"),width = 4,height = 3)}
}



###highlight Cluster13
HSMM@colData$Cluster13_status<-"Neg"
HSMM@colData$Cluster13_status[HSMM@colData$refinedID=="Cluster_3"]="Pos"

plot_cells(HSMM,color_cells_by = "Cluster13_status", group_label_size = 3,cell_size = 1,show_trajectory_graph = F)



