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
setwd("/Users/hsong/Desktop/Bladder_2023/Epithelial/")

###Load in the final object
dge<-readRDS("../Bladder_all_Final.rds")
cells_endo<-colnames(dge)[dge$Final_ID=="Endothelial"]
cells_bcell<-colnames(dge)[dge$Final_ID=="Bcell"]
cells_fibro<-colnames(dge)[dge$Final_ID=="Fibroblast"]
cells_sm<-colnames(dge)[dge$Final_ID=="Smooth_Muscle"]
cells_mast<-colnames(dge)[dge$Final_ID=="Mast_Cell"]
cells_myeloid<-colnames(dge)[dge$Final_ID=="Myeloid"]
cells_neutrophil<-colnames(dge)[dge$Final_ID=="Neutrophil"]
cells_plasma<-colnames(dge)[dge$Final_ID=="Plasma"]
cells_tcell<-colnames(dge)[dge$Final_ID=="T-cell"]

cells_epithelial<-colnames(Epithelial)
### subset only the epithelial cells
dge <- subset(dge,cells = cells_epithelial)
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


####Sankey diagram
library(networkD3)

nodes=data.frame("name" = c(unique(dge$Name2),unique(dge$Name)))

dfr<-as.data.frame(table(dge$Name2,dge$Name))
names(dfr)<-c("source","target","value")

## create a dataframe with 10 nodes
nodes=unique(c(as.character(dfr$source), as.character(dfr$target)))
dfr_link <- data.frame(source=as.numeric(factor(dfr$source, levels=nodes))-1, target=as.numeric(factor(dfr$target, levels=nodes))-1, value=dfr$value)
nodes=as.data.frame(nodes)
colnames(nodes)="name"
p = sankeyNetwork(Links = dfr_link, Nodes = nodes,
                  Source = "source", Target = "target",
                  Value = "value", NodeID = "name",
                  fontSize = 16, nodeWidth = 40)
p






