
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

####Dividing cells
dge<-Tumor_variant
dge$Dividing<-"Non_Dividing"
dge$Dividing[dge@active.ident %in% c("MicropapillaryJM_Cluster_1","MicropapillaryPS_Cluster_3",
                                     "Micropapillary21222_Cluster_2","Plasmacytoid11734_Cluster_2",
                                     "Small_Large_Cluster_1","Small_Large_Cluster_3",
                                     "Squamousdiff21226_Cluster_2","PureUC12923_Cluster_4",
                                     "PureUCSG_Cluster_4","PureUCHG_Cluster_4")]="Dividing"
dge<-SetIdent(dge,value = as.vector(dge$Dividing))
Volcano_Plot(dge,"Dividing_NonDividing",30)
Gene_Dividing<-c("TOP2A","CENPF","MKI67","HMGB2","TPX2",
                 "ASPM","KIF22","PRC1","ATAD2","SMC4",
                 "NUF2","ANLN","CENPE","STMN1","BIRC5",
                 "NCAPG","CDK1","KPNA2","H2AFZ","NUSAP1",
                 "TUBB2B","CASC5","TUBA1A","KIF11","KIF23")

IDs_Dividing<-c("MicropapillaryJM_Cluster_1","MicropapillaryPS_Cluster_3",
                "Micropapillary21222_Cluster_2","Plasmacytoid11734_Cluster_2",
                "Small_Large_Cluster_1","Small_Large_Cluster_3",
                "Squamousdiff21226_Cluster_2","PureUC12923_Cluster_4",
                "PureUCSG_Cluster_4","PureUCHG_Cluster_4")

dge_temp<-subset(dge_temp,downsample=200)

dge<-SetIdent(dge,value = as.vector(dge$ID_subcluster))
sub_ave<-AverageExpression(dge,return.seurat = T)
dge_temp<-sub_ave
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
z <- DoHeatmap(dge_temp, features = Gene_Dividing) + NoLegend()
y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()

# my_sample_col <- data.frame((dge_temp@active.ident))
#my_sample_col <- data.frame((dge$cell_type_annotation))
#row.names(my_sample_col) <- colnames(data_subset)
dge_temp$Dividing="Non_Dividing"
dge_temp$Dividing[colnames(dge_temp) %in% IDs_Dividing]="Dividing"
my_sample_col <- data.frame((dge_temp$Dividing))
row.names(my_sample_col) <- colnames(dge_temp)
colnames(my_sample_col)="Annotation"

ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = T,clustering_method="ward.D2",cluster_cols = T,cluster_rows = F)
save_plot_pdf(ww, "Heatmap_clustered_nopureUC_Para.pdf",15,10)

w<-w[,order(colnames(w))]
w<-w[,c(32:46,6:31,1:5,47:55)]
ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = T,clustering_method="ward.D2",cluster_cols = F,cluster_rows = F)
save_plot_pdf(ww, "Heatmap_unclustered_nopureUC_Para.pdf",15,10)

ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = T,clustering_method="ward.D2",cluster_cols = T,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_clustered_nopureUC_Para.pdf",15,10)

w<-w[,order(colnames(w))]
w<-w[,c(29,34:48,6:28,31:33,1:5,49:57,30)]
ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = T,clustering_method="ward.D2",cluster_cols = F,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_unclustered_nopureUC_Para.pdf",15,10)






