library("Hmisc")
library(corrplot)
#### Identify correlated genes with TM4SF1 in Tumor cells
options(warn = 1)
Epithelial<-readRDS("../Epithelial.rds")
ID="Variant"
dge<-subset(Epithelial,idents=c("Variant"))
dge<-subset(Epithelial,idents=c("Pure_UC"))
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- ScaleData(object=dge,features=rownames(dge))

matrix<-dge@assays$RNA@data
matrix_mod<-as.matrix(matrix)

gene<-as.numeric(matrix_mod["TM4SF1",])
correlations<-apply(matrix_mod,1,function(x){rcorr(gene,x)})
DF=NULL
for (i in 1:length(correlations)){
  DF$Gene[i]=names(correlations)[i]
  DF$Correlation[i]=correlations[[i]]$r[1,2]
  DF$P_value[i]=correlations[[i]]$P[1,2]
}
DF<-as.data.frame(DF)
#View(DF)

correlations<-DF
correlations$Cluster=ID
correlations<-correlations[order(correlations$P_value,decreasing = F),]
write.table(correlations,paste0("Correlation_TM4SF1_",ID,".txt"),sep="\t",col.names = T,row.names = T,quote = F)


Tumor_DF<-read.table("./Correlation_TM4SF1_Variant.txt",sep="\t",header = T)
Tumor_DF<-Tumor_DF[Tumor_DF$P_value<0.05,]
Tumor_DF<-Tumor_DF[abs(Tumor_DF$Correlation)>=0.2,]
Tumor_DF<-Tumor_DF[!is.na(Tumor_DF$P_value),]
Tumor_DF<-Tumor_DF[order(Tumor_DF$Correlation,decreasing = T),]
Tumor_DF<-Tumor_DF[-c(1),]

NonTumor_DF<-read.table("./Correlation_TM4SF1_Pure UC.txt",sep="\t",header = T)


DF<-Tumor_DF
DF$Correlation_UC=0
DF$Pvalue_UC=1
for (i in 1:nrow(DF)){

DF$Correlation_UC[i]<-NonTumor_DF$Correlation[NonTumor_DF$Gene==DF$Gene[i]]
DF$Pvalue_UC[i]<-NonTumor_DF$P_value[NonTumor_DF$Gene==DF$Gene[i]]
}
colnames(DF)[2]="Correlation_Variant"
colnames(DF)[3]="Pvalue_Variant"
DF<-DF[,-c(4)]


matrix_mod<-as.matrix(matrix)
matrix_mod<-matrix_mod[(rownames(matrix_mod) %in% NonTumor_DF$Gene),]

c_df <- rcorr(as.matrix(t(matrix_mod)))
View(c_df)

pdf("Correlation_NonTumor.pdf")
corrplot(corr=c_df$r, p.mat=c_df$P, sig.level=0.05, 
         method='color', diag=FALSE, addCoef.col=1, type='upper', insig='blank',
         number.cex=.8)
dev.off()

TM4SF1_Correlation_Tumoronly<-Tumor_DF[!(Tumor_DF$Gene %in% NonTumor_DF$Gene),]
write.table(TM4SF1_Correlation_Tumoronly,"TM4SF1_Correlation_Variantonly.txt",sep="\t",col.names = T,row.names = F,quote = F)

TM4SF1_Correlation_Tumoronly<-read.table("./TM4SF1_Correlation_Variantonly.txt",sep = "\t",header = T)
sig_cluster<-TM4SF1_Correlation_Tumoronly$Gene
###Make heatmap with dendrogram

library(tidyr)
library(textshape)
genelist<-sig_cluster

dge_temp<-SetIdent(Epithelial,value = as.vector(Epithelial$histology2))

dge_temp<-subset(dge_temp,cells=colnames(dge_temp)[!(dge_temp$histology2 %in% c("Paired Normal","Paraganglioma"))])
Celltype<-unique(dge_temp$histology2)
Celltype<-Celltype[order(Celltype)]
Celltype<-Celltype[c(10,9,8,1:7,11:12)]
levels(dge_temp@active.ident)=Celltype
dge_temp<-subset(dge_temp,downsample=200)
dge_temp<-AverageExpression(dge,return.seurat = T)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
genelist<-c("PROM1","POU5F1","SOX2", "ALDH1A1", "SOX4", "EZH2", "YAP1", "CD44", "KRT14")
z <- DoHeatmap(dge_temp, features = genelist) + NoLegend()
y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()

my_sample_col <- data.frame((dge_temp@active.ident))

colnames(my_sample_col)="Annotation"


ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = T,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_clustered_nopureUC_Para.pdf",20,30)

ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = F,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_unclustered_nopureUC_Para.pdf",20,30)

ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = T,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_clustered_detail_nopureUC_Para.pdf",20,30)

ww<-pheatmap::pheatmap(w,annotation_col = my_sample_col,show_colnames = F,clustering_method="ward.D2",cluster_cols = F,cluster_rows = T)
save_plot_pdf(ww, "Heatmap_unclustered_detail_nopureUC_Para.pdf",20,30)


####Correlation Volcano plot
source("~/Desktop/Volcano_Plot.R")
dge$Correlation_Cluster="Upper"
dge$Correlation_Cluster[dge@active.ident %in% c("PureUCHG_Cluster_2","PureUC12923_Cluster_1","PureUC12923_Cluster_4",
                                                "Small_Large_Cluster_4","Small_Large_Cluster_1","Small_Large_Cluster_3",
                                                "Squamousdiff21226_Cluster_3","Lymphoepithelioma20847_Cluster_1","Squamousdiff21226_Cluster_1",
                                                "Squamousdiff21226_Cluster_2","Nested12041_Cluster_2","Nested12041_Cluster_0",
                                                "Nested12041_Cluster_1","Nested12041_Cluster_3",
                                                "Plasmacytoid11734_Cluster_0","Plasmacytoid11734_Cluster_1","Plasmacytoid11734_Cluster_2")]="Lower"
dge<-SetIdent(dge,value = as.vector(dge$Correlation_Cluster))
Volcano_Plot(dge,"Correlation_Heatmap_Volcano",30,20,20)


####GO analysis of Correlation_based clusters
source("~/Desktop/GO_analysis.R")
GO_analysis(dge,"Correlation_based",0.5)
#### Load in the Signatuer gene sets and confirm the signature scores in the violin plots
setwd("./Correlation_GO_analysis/")
List_cluster<-read.table("./List_Cluster.txt",sep="\n",header = F)[[1]]
dge<-SetIdent(dge,value = as.vector(dge$ID_subcluster))
for (i in 1:length(List_cluster)){
  Signature <- read.table(paste0("Correlation_based_",List_cluster[i],"_Signature.txt"),sep="\n",header = F)[[1]]
  Name <- List_cluster[i]
  dge <- AddModuleScore(object = dge,assay = "RNA", features = list(Signature), ctrl = 5, name = Name)
  V1<-as.vector(FetchData(object = dge,vars = paste0(Name,"1")))[[1]]
  BOX_df<-NULL
  BOX_df$id<-as.vector(dge@active.ident)
  BOX_df$value<-as.vector(V1)
  BOX_df<-data.frame(BOX_df)
  FeaturePlot(dge,features=paste0(c(Name),"1"))+scale_color_gradientn( colours = c("grey","red"),  limits = c(0, max(V1)))
  ggsave(file=paste0("Modulescore_",Name,"_feature.pdf"),width = 15,height = 15,units = "cm")
  
  
  ggplot(BOX_df, aes(x=reorder(id,-value), y=value,fill=id)) + geom_boxplot() + #geom_jitter(shape=16,position = position_jitter(0.1))+
    #stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ 
    theme(legend.position="none",text = element_text(size=6))+
    stat_compare_means(method = "anova",label.x = 3,label.y = max(unlist(V1))+0.05)+
    ggtitle(Name)+rotate_x_text(angle = 45)
  ggsave(file=paste0("Modulescore_",Name,"_Boxplot.pdf"),width = 40,height = 15,units = "cm")
  
}

