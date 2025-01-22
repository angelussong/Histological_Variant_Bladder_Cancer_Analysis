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

sce <- SingleCellExperiment(list(counts = dge@assays$RNA@counts))
sce <- decontX(sce)
sce_count<-decontXcounts(sce)
dge<- CreateSeuratObject(counts = sce_count,project = "Bladder_Final", min.cells = 0,min.features = 0,meta.data = dge@meta.data)

dge$Name="Unknown"
dge$Name[dge$patient %in% c("21217_N", "21217_T")]="VAR11"
dge$Name[dge$patient %in% c("21222_N", "21222_T")]="VAR1"
dge$Name[dge$patient %in% c("21226_N", "21226_T")]="VAR7"
dge$Name[dge$patient %in% c("21262_T")]="VAR9"
dge$Name[dge$patient %in% c("FG_CIS","FG_N","FG_T")]="VAR5"
dge$Name[dge$patient %in% c("HG_T")]="UC3"
dge$Name[dge$patient %in% c("20847_T")]="VAR6"
dge$Name[dge$patient %in% c("21032_N","21032_T")]="VAR12"
dge$Name[dge$patient %in% c("JM_T","JMx_N")]="VAR3"
dge$Name[dge$patient %in% c("PG_SN","PG_MT","PG_ST")]="PG"
dge$Name[dge$patient %in% c("PS_T")]="VAR2"
dge$Name[dge$patient %in% c("11734")]="VAR8"
dge$Name[dge$patient %in% c("12041")]="VAR4"
dge$Name[dge$patient %in% c("12050_N","12050_T")]="VAR10"
dge$Name[dge$patient %in% c("12923")]="UC1"
dge$Name[dge$patient %in% c("SG")]="UC2"

dge<-subset(dge,cells=colnames(dge)[!(dge$patient=="RC_N")])
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- FindVariableFeatures(object = dge, selection.method = "vst", nfeatures =3000)

# # Identify the 10 most highly variable genes
top10 <- head(x = VariableFeatures(object = dge), 10)
top10
# 
# # plot variable features with and without labels
plot1 <- VariableFeaturePlot(object = dge)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE,xnudge = 0,ynudge = 0)
CombinePlots(plots = list(plot1, plot2))
ggsave(file="VariableFeature_3000.pdf",width = 20,height = 20,units = "cm")
# 
# all.genes <- rownames(x = dge)
dge <- ScaleData(object=dge,features=VariableFeatures(object = dge))
dge <- RunPCA(object = dge, features = VariableFeatures(object=dge),npcs = 100)


# Create a vector of gene names for downstream analysis
#gene.names.vector <- as.character("gene names")
# Run PCA on selected genes
#dge <- RunPCA(object = dge, features = gene.names.vector, do.print = TRUE, pcs.print = 1:5, 
#                 genes.print = 5)

slot(dge[["pca"]], "misc")
print(x = dge[["pca"]], dims = 1:5, nfeatures = 5)
mat <- Seurat::GetAssayData(dge,assay="RNA",slot="scale.data")
pca <-dge[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))
eigValues = (pca@stdev)^2
varExplained = eigValues / total_variance
sum(varExplained)

DimPlot(object = dge, reduction = "pca",group.by = "ID2")
ggsave(file="PCA_byID.pdf",width = 25,height = 20,units = "cm")
DimPlot(object = dge, reduction = "pca",group.by = "Name")
ggsave(file="PCA_byName.pdf",width = 20,height = 20,units = "cm")

png("PCA_heatmap.png")
DimHeatmap(object = dge, dims = 1:10, cells = length(dge$nCount_RNA), balanced = TRUE)
dev.off()
png("Elbowplot_dge.png")
ElbowPlot(dge,ndims = 100)
dev.off()

n_pc=100
dge <- FindNeighbors(dge, dims = 1:n_pc,k.param = 30)
dge <- FindClusters(dge, resolution =c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

library(clustree)
###epithelial
pdf("Stability.pdf",width = 15,height = 25)
clustree(dge, prefix = "RNA_snn_res.")
dev.off()

dge <- FindClusters(dge, resolution =0.4)

dge=RunTSNE(dge,dims.use = 1:n_pc,perplexity=40,seed.use = 10,check_duplicates = F)
TSNEPlot(dge,pt.size = 1,label=T)+NoLegend()
ggsave(file="TSNE_raw.pdf",width = 20,height = 20,units = "cm")

dge <- RunUMAP(dge, dims = 1:n_pc)
DimPlot(dge, reduction = "umap",label=T)+NoLegend()
ggsave(file="Umap_raw.pdf",width = 20,height = 20,units = "cm")


DimPlot(dge, reduction = "umap",group.by="Name",label=F)
ggsave(file="Umap_bypatient.pdf",width = 25,height = 20,units = "cm")

# Print number of cells per cluster
Tab<-table(dge$orig.ident,Idents(object = dge))
write.table(Tab,"table_initialclustering.txt",sep="\t",row.names = T,col.names = T)
#clustering dendrogram
dge<- BuildClusterTree(dge)
png("Cluster_tree_integrated.png")
PlotClusterTree(dge)
dev.off()

dge$ID<-dge$predicted.id
dge<-SetIdent(dge,value = as.vector(dge$predicted.id))
dge_temp<-subset(dge,downsample=200)
pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
write.table(pbmc.markers,"allcell_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
write.table(pbmc.markers,"Tumor_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
write.table(pbmc.markers,"allcell_markers_byCellType.txt",col.names = T,row.names = T,quote = F,sep="\t")

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top5,"top5_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
write.table(top10,"top10_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")
write.table(top10,"top10_markers_bytype.txt",col.names = T,row.names = T,quote = F,sep="\t")

dge_temp<-subset(dge,downsample=200)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp, features = top10$gene,raster = F) + theme(axis.text.y = element_text(size = 5))+NoLegend()
ggsave(file="Heatmap_10_bycluster.png",width = 30,height = 30,units = "cm",limitsize = FALSE)
ggsave(file="Heatmap_10_bycelltype.png",width = 30,height = 30,units = "cm",limitsize = FALSE)

dge <- ScaleData(object=dge,features=rownames(dge))
DoHeatmap(dge, features = top10$gene,raster = F,size = 2) + theme(axis.text.y = element_text(size = 5))+NoLegend()
ggsave(file="Heatmap_10_Tumor.png",width = 15,height = 15,units = "cm",limitsize = FALSE)


DotPlot(dge, features = rev(unique(top10$gene)), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave(file="Dotplot_10_bytype.pdf",width = 16,height = 6)


###Transfer labels
transfer.anchors <- FindTransferAnchors(
  reference = dge_ref,
  query = dge,
  reduction = 'pcaproject'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = dge_ref$toptier_annotation,
  # weight.reduction = dge[['lsi']],
  # dims = 2:50
)

dge <- AddMetaData(object = dge, metadata = predicted.labels)
cells_tmp<-colnames(dge)[(dge$seurat_clusters=="11") & (dge$predicted.id=="Smooth_Muscle")]
cells_tmp2<-colnames(dge)[dge@reductions$umap@cell.embeddings[,1]>0]
cells_tmp<-intersect(cells_tmp,cells_tmp2)
dge$Final_ID<-dge$predicted.id
dge$Final_ID[dge$seurat_clusters %in% c("16","17","20","31")]="Endothelial"
dge$Final_ID[colnames(dge) %in% cells_tmp]="Endothelial"
cells_tmp<-colnames(dge)[(dge$seurat_clusters=="11") & (dge$predicted.id=="Fibroblast")]
cells_tmp2<-colnames(dge)[dge@reductions$umap@cell.embeddings[,1]<0]
cells_tmp<-intersect(cells_tmp,cells_tmp2)
dge$Final_ID[colnames(dge) %in% cells_tmp]="Smooth_Muscle"
dge$Final_ID[dge$seurat_clusters %in% c("14","18","25","33","30")]="Epithelial"
dge$Final_ID[dge$seurat_clusters %in% c("3")]="T-cell"
dge$Final_ID[dge$seurat_clusters %in% c("4","6","9","24","27")]="Myeloid"
dge$Final_ID[dge$seurat_clusters %in% c("29")]="Neutrophil"
dge$Final_ID[dge$seurat_clusters %in% c("22")]="Mast_Cell"

dge<-SetIdent(dge,value = as.vector(dge$Final_ID))


###prepare invferCNV
dge$InferCNV_ID<-as.character(dge$seurat_clusters)
dge$InferCNV_ID[dge$Final_ID=="Endothelial"]="Endothelial"
dge$InferCNV_ID[dge$Final_ID=="Fibroblast"]="Fibroblast"
dge$InferCNV_ID[dge$Final_ID=="Mast_Cell"]="Mast_Cell"
dge$InferCNV_ID[dge$Final_ID=="Myeloid"]="Myeloid"
dge$InferCNV_ID[dge$Final_ID=="Neutrophil"]="Neutrophil"
dge$InferCNV_ID[dge$Final_ID=="Smooth_Muscle"]="Smooth_Muscle"
dge$InferCNV_ID[dge$Final_ID=="T-cell"]="T-cell"
dge$InferCNV_ID[dge$InferCNV_ID=="11"]="Smooth_Muscle"

dge$Name[dge$Name=="VAR1"]="VAR01"
dge$Name[dge$Name=="VAR2"]="VAR02"
dge$Name[dge$Name=="VAR3"]="VAR03"
dge$Name[dge$Name=="VAR4"]="VAR04"
dge$Name[dge$Name=="VAR5"]="VAR05"
dge$Name[dge$Name=="VAR6"]="VAR06"
dge$Name[dge$Name=="VAR7"]="VAR07"
dge$Name[dge$Name=="VAR8"]="VAR08"
dge$Name[dge$Name=="VAR9"]="VAR09"
dge$Name[dge$Name=="UC1"]="UC01"
dge$Name[dge$Name=="UC2"]="UC02"
dge$Name[dge$Name=="UC3"]="UC03"
dge$Name[dge$Name=="VAR12"]="UC04"
dge$InferCNV_ID<-dge$Name
dge$InferCNV_ID[dge$Final_ID=="Endothelial"]="Endothelial"
dge$InferCNV_ID[dge$Final_ID=="Fibroblast"]="Fibroblast"
dge$InferCNV_ID[dge$Final_ID=="Mast_Cell"]="Mast_Cell"
dge$InferCNV_ID[dge$Final_ID=="Myeloid"]="Myeloid"
dge$InferCNV_ID[dge$Final_ID=="Neutrophil"]="Neutrophil"
dge$InferCNV_ID[dge$Final_ID=="Smooth_Muscle"]="Smooth_Muscle"
dge$InferCNV_ID[dge$Final_ID=="T-cell"]="T-cell"

dge<-SetIdent(dge,value = as.vector(dge$InferCNV_ID))
dge<-SetIdent(dge,value = as.vector(dge$InferCNV_ID2))
dge_temp<-subset(dge,downsample=500)
dge_infer_mtx <- as.matrix(GetAssayData(dge_temp, slot = "counts"))
write.table(dge_infer_mtx,"Bladder_patient_counts.txt",sep="\t",row.names = T,col.names = T,quote=F)
phenotype<-NULL
phenotype$cellname<-colnames(dge_temp)
phenotype$ID<-as.vector(dge_temp@active.ident)
phenotype<-as.data.frame(phenotype)
#View(phenotype)
write.table(phenotype,"Bladder_patient_phenotype.txt",sep="\t",row.names = F,col.names = F,quote=F)

Featurename<-c("BE","Endothelial","Fib","SM","Tcell","Myeloid","Neutrophil","Mast_cell")
Features<-list()
Features[[1]] <- read.table("~/Desktop/Seqwell_combined/TCGA/BE_Markers.txt",sep="\n",header = F)[[1]]
Features[[2]]<-read.table("~/Desktop/Seqwell_combined/TCGA/Endothelial_Markers.txt",sep="\n",header = F)[[1]]
Features[[3]]<-read.table("~/Desktop/Seqwell_combined/TCGA/Fibroblast_Markers.txt",sep="\n",header = F)[[1]]
Features[[4]]<-read.table("~/Desktop/Seqwell_combined/TCGA/Smooth_muscle_Markers.txt",sep="\n",header = F)[[1]]
Features[[5]]<-read.table("~/Desktop/Seqwell_combined/TCGA/T-cell_Markers.txt",sep="\n",header = F)[[1]]
Features[[6]]<-read.table("~/Desktop/Seqwell_combined/TCGA/Myeloid_Markers.txt",sep="\n",header = F)[[1]]
Features[[7]]<-c("FCGR3B","IL8","FPR1","TNFAIP6","TREM1")
Features[[8]]<-c("CPA3","KIT","MS4A2","TPSAB1")

for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  dge <- AddModuleScore(object = dge, features = list(gene.set), ctrl = 2, name = Featurename[i])
  V1<-as.vector(FetchData(object = dge,vars = paste0(Featurename[i],"1")))
  FeaturePlot(dge,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(unlist(V1))))
  ggsave(file=paste0("Featureplot_",Featurename[i],".pdf"),width = 8,height = 8)
  BOX_df<-NULL
  BOX_df$id<-dge@active.ident
  BOX_df$value<-as.vector(V1)
  BOX_df<-data.frame(BOX_df)
  colnames(BOX_df)<-c("id","value")
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot() + #geom_violin() +#geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=3,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
    stat_compare_means(method = "anova",label.x = 3,label.y = max(unlist(V1))+0.05)+
    ggtitle(Featurename[i])+rotate_x_text(angle = 45)
  ggsave(file=paste0("Vlnplot_",Featurename[i],".pdf"),width = 10,height = 5)
  
}

VlnPlot(dge,features=c("PIGR","LTF","MMP7","SCGB1A1"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="Club_markers_bycluster.png",width = 35,height = 20,units = "cm")

VlnPlot(dge,features=c("AMACR","PCA3","ERG","ETV1"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="Tumor_markers_bycluster.png",width = 35,height = 20,units = "cm")
VlnPlot(dge,features=c("PCAT4","PCAT14","SPON2","GOLM1"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="Tumor_markers2_bycluster.png",width = 35,height = 20,units = "cm")

VlnPlot(dge,features=c("KLK3","KRT8","AR","NKX3-1"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="LE_markers_bycluster.png",width = 35,height = 20,units = "cm")

VlnPlot(dge,features=c("KRT5","KRT15","DST","TP63"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="BE_markers_bycluster.png",width = 35,height = 20,units = "cm")

VlnPlot(dge,features=c("EPCAM","PTPRC","IL7R","MKI67"),ncol=2)+NoLegend()+rotate_x_text(angle = 45)
ggsave(file="Misc_markers_bycluster.png",width = 35,height = 20,units = "cm")

FeaturePlot(dge,features=c("PIGR","LTF","MMP7","SCGB1A1"),ncol=2)
ggsave(file="Club_markers_feature.png",width = 30,height = 30,units = "cm")

FeaturePlot(dge,features=c("AMACR","PCA3","ERG","ETV1"),ncol=2)
ggsave(file="Tumor_markers_feature.png",width = 30,height = 30,units = "cm")

FeaturePlot(dge,features=c("KLK3","KRT8","AR","NKX3-1"),ncol=2)
ggsave(file="LE_markers_feature.png",width = 30,height = 30,units = "cm")

FeaturePlot(dge,features=c("KRT5","KRT15","DST","TP63"),ncol=2)
ggsave(file="BE_markers_feature.png",width = 30,height = 30,units = "cm")

FeaturePlot(dge,features=c("EPCAM","PTPRC","IL7R","MKI67"),ncol=2)
ggsave(file="Misc_markers_feature.png",width = 30,height = 30,units = "cm")



###Remove paired normal
dge$Final_ID[dge$seurat_clusters %in% c("18")]="Plasma"
dge$Final_ID[dge$seurat_clusters %in% c("25")]="Bcell"
dge$Final_ID[dge$seurat_clusters %in% c("33")]="Myeloid"
unique(dge$patient)
dge<-subset(dge,cells=colnames(dge)[!dge$patient %in% c("21217_N","21222_N","21226_N","FG_N","21032_N","JMx_N","PG_SN","12050_N")])
table(dge$Name,dge$Final_ID)
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

unique(dge$toptier_annotation)
DimPlot(dge,group.by = "toptier_annotation")


##
mtx<-read.table("./coadread_tcga_pan_can_atlas_2018/data_mrna_seq_v2_rsem.txt")
