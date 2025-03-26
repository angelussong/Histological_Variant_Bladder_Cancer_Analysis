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

####Box plots of multiple genesets
KEGG_CAM1<-read.table("../Draft/KEGG_CAMs.csv",sep="\n",header = F)[[1]]
KEGG_CAM2<-read.table("../Draft/KEGG_CAMs1.csv",sep="\n",header = F)[[1]]
Leu_Adhesion1<-read.table("../Draft/Leukocyte_adhesion.csv",sep="\n",header = F)[[1]]
Leu_Adhesion2<-read.table("../Draft/Leukocyte_adhesion1.csv",sep="\n",header = F)[[1]]


Name="Leu_Adhesion2"
gene.set <- Leu_Adhesion2
dge <- AddModuleScore(object = dge,assay = "RNA", features = list(gene.set), ctrl = 5, name = Name)
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

####Boxplot in the TCGA
TCGA<-readRDS("~/Desktop/Bladder/blca_tcga_pan_can_atlas_2018/TCGA_data.rds")
Name="Leu_Adhesion2"
gene.set <- Leu_Adhesion2
TCGA <- AddModuleScore(object = TCGA,assay = "RNA", features = list(gene.set), ctrl = 5, name = Name)
V1<-as.vector(FetchData(object = TCGA,vars = paste0(Name,"1")))[[1]]
V1<-as.vector(FetchData(object = TCGA,vars = "NRXN3"))[[1]]
BOX_df<-NULL
BOX_df$id<-as.vector(TCGA@active.ident)
BOX_df$value<-as.vector(V1)
BOX_df<-data.frame(BOX_df)

ggplot(BOX_df, aes(x=reorder(id,-value), y=value,fill=id)) + geom_boxplot() + geom_jitter(shape=16,position = position_jitter(0.1))+
  #stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+ 
  theme(legend.position="none",text = element_text(size=6))+
  stat_compare_means(method = "anova",label.x = 3,label.y = max(unlist(V1))+0.05)+
  ggtitle(Name)+rotate_x_text(angle = 45)
ggsave(file=paste0("TCGA_Modulescore_",Name,"_Boxplot.pdf"),width = 20,height = 15,units = "cm")


####Derive Specific CAM genes
dge<-readRDS("./dge_epithelial3.rds")
dge<-SetIdent(dge,value = as.vector(dge$histology2))
Gene_list<-read.table("./KEGG_CAMs.csv",sep="\n",header = F)[[1]]
Gene_list<-intersect(Gene_list,rownames(dge))
Cluster_ID<-as.vector(unique(dge@active.ident))
for (i in 1:length(Cluster_ID)){
  DF<-NULL
  DF$Gene<-Gene_list
  DF$log2FC=0
  DF$p_value=1
  DF$FDR=1
  DF$Percentage=0
  DF$Expression=0
  DF<-as.data.frame(DF)
  dge_target<-subset(dge,idents=Cluster_ID[i])
  dge_rest<-subset(dge,idents=Cluster_ID[i],invert=T)
  for (j in 1:length(Gene_list)){
    Target_expression<-as.numeric(FetchData(dge_target,vars = c(Gene_list[j]))[,1])
    Rest_expression<-as.numeric(FetchData(dge_rest,vars = c(Gene_list[j]))[,1])
    
    DF$log2FC[j]=log2(mean(Target_expression)/mean(Rest_expression))
    wilcox_test<-wilcox.test(Target_expression,Rest_expression,alternative = c("two.sided"),paired = F)
    DF$p_value[j]=wilcox_test$p.value
    DF$Expression[j]=mean(Target_expression)
    if(mean(Target_expression)>0){
      DF$Percentage[j]=length(Target_expression[Target_expression>0])/ncol(dge_target)*100
    }
  }
  DF$FDR<-p.adjust(DF$p_value,n = nrow(DF),method = "fdr")
  DF<-DF[order(DF$log2FC,decreasing = T),]
  DF<-DF[DF$FDR<0.05,]
  DF<-DF[DF$log2FC>=1,]
  DF<-DF[!is.na(DF$Gene),]
  write.table(DF,paste0(Cluster_ID[i],"_KEGG_CAM_Genelist.txt"),sep="\t",row.names = F,col.names = T,quote = F)
}



####Correlation Heatmap Extension
dge<-readRDS("./dge_epithelial3.rds")

dge$Correlation_cluster="Group1"
dge$Correlation_cluster[dge$ID_subcluster %in% c("PureUCSG_2")]="Group2_PureUCSG_2"
dge$Correlation_cluster[dge$ID_subcluster %in% c("PureUC12923_1","PureUC12923_4")]="Group2_PureUC12923_1_4"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Small_Large_1","Small_Large_3","Small_Large_4")]="Group2_Small_Large"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Squamousdiff21226_3")]="Group2_Squamous21226_3"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Lymphoepithelioma20847_1")]="Group2_Lympho20847_1"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Squamousdiff21226_1","Squamousdiff21226_2")]="Group2_Squamous21226_1_2"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Nested12041_1","Nested12041_2","Nested12041_3","Nested12041_0")]="Group2_Nested12041"
dge$Correlation_cluster[dge$ID_subcluster %in% c("Plasmacytoid11734_2","Plasmacytoid11734_0","Plasmacytoid11734_1")]="Group2_Plasmacytoid11734"


dge<-SetIdent(dge,value = as.vector(dge$Correlation_cluster))
pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes,RPL.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
write.table(pbmc.markers,"Correlation_based_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.table(top10,"top10_Correlation_markers.txt",col.names = T,row.names = T,quote = F,sep="\t")

dge_temp<-subset(dge,downsample=100)
dge_temp <- ScaleData(object=dge_temp,features=rownames(dge_temp))
DoHeatmap(dge_temp, features = top10$gene,raster = F) + theme(axis.text.y = element_text(size = 5))+NoLegend()
ggsave(file="Heatmap_10.png",width = 25,height = 30,units = "cm",limitsize = FALSE)
# can use a RColorBrewer pallete, I chose RdBu here
DotPlot(dge, features = rev(unique(top10$gene)), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))

ggsave(file="dot_plot_Correlation_based.pdf",width = 45,height = 15,units = "cm",limitsize = FALSE)

