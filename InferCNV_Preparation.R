library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(clustree)
library(celda)
library(SingleCellExperiment)
library(jjb)

theme_set(theme_cowplot())
setwd("/Users/hsong/Desktop/Bladder_2023/")
source("~/Desktop/save_pheatmap_pdf.R")
source("~/Desktop/save_plot_pdf.R")
source("~/Desktop/Convert_Seurat_Scanpy.R")
source("~/Desktop/Volcano_Plot.R")
source("~/Desktop/GO_analysis.R")
source("~/Desktop/Run_Fisher.R")
source("~/Desktop/Convert_rawcount.R")

dge<-readRDS("./dge_mergedall_032123.rds")

###prepare invferCNV
dge$InferCNV_ID2<-as.character(dge$seurat_clusters)
dge$InferCNV_ID2[dge$Final_ID=="Endothelial"]="Endothelial"
dge$InferCNV_ID2[dge$Final_ID=="Fibroblast"]="Fibroblast"
dge$InferCNV_ID2[dge$Final_ID=="Mast_Cell"]="Mast_Cell"
dge$InferCNV_ID2[dge$Final_ID=="Myeloid"]="Myeloid"
dge$InferCNV_ID2[dge$Final_ID=="Neutrophil"]="Neutrophil"
dge$InferCNV_ID2[dge$Final_ID=="Smooth_Muscle"]="Smooth_Muscle"
dge$InferCNV_ID2[dge$Final_ID=="T-cell"]="T-cell"
dge$InferCNV_ID2[dge$InferCNV_ID=="11"]="Smooth_Muscle"

dge$InferCNV_ID<-dge$Name
dge$InferCNV_ID[dge$Final_ID=="Endothelial"]="Endothelial"
dge$InferCNV_ID[dge$Final_ID=="Fibroblast"]="Fibroblast"
dge$InferCNV_ID[dge$Final_ID=="Mast_Cell"]="Mast_Cell"
dge$InferCNV_ID[dge$Final_ID=="Myeloid"]="Myeloid"
dge$InferCNV_ID[dge$Final_ID=="Neutrophil"]="Neutrophil"
dge$InferCNV_ID[dge$Final_ID=="Smooth_Muscle"]="Smooth_Muscle"
dge$InferCNV_ID[dge$Final_ID=="T-cell"]="T-cell"
dge$InferCNV_ID[dge$Final_ID=="Plasma"]="Plasma"
dge$InferCNV_ID[dge$Final_ID=="Bcell"]="Bcell"


dge<-SetIdent(dge,value = as.vector(dge$InferCNV_ID))
dge<-SetIdent(dge,value = as.vector(dge$InferCNV_ID2))
dge_temp<-subset(dge,cells=colnames(dge)[!(dge$InferCNV_ID %in% c("Mast_Cell","Neutrophil","Bcell"))])
dge_temp<-SetIdent(dge_temp,value = as.vector(dge_temp$InferCNV_ID))

dge_temp<-subset(dge,downsample=500)
dge_infer_mtx <- as.matrix(GetAssayData(dge_temp, slot = "counts"))
write.table(dge_infer_mtx,"Bladder_patient_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)
phenotype<-NULL
phenotype$cellname<-colnames(dge_temp)
phenotype$ID<-as.vector(dge_temp@active.ident)
phenotype<-as.data.frame(phenotype)
#View(phenotype)
write.table(phenotype,"Bladder_patient_phenotype.txt",sep="\t",row.names = F,col.names = F,quote=F)


