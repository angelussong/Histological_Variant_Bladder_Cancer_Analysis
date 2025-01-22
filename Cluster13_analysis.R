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
setwd("/Users/hsong/Desktop/Bladder_2023/Cluster13/")
source("~/Desktop/save_pheatmap_pdf.R")
source("~/Desktop/save_plot_pdf.R")
source("~/Desktop/Convert_Seurat_Scanpy.R")
source("~/Desktop/Volcano_Plot.R")
source("~/Desktop/GO_analysis.R")
source("~/Desktop/Run_Fisher.R")
source("~/Desktop/Convert_rawcount.R")

dge<-readRDS("./dge_13.rds")
Cluster13_markers<-read.table("~/Desktop/Bladder_2023/Cluster13/epithelial_allcell_markers_bycluster.txt",sep="\t",header = T)
Cluster13_markers<-Cluster13_markers[Cluster13_markers$cluster=="Cluster13",]
Cluster13_markers<-Cluster13_markers[Cluster13_markers$avg_log2FC>1,]
Cluster13_markers<-Cluster13_markers[Cluster13_markers$p_val_adj<0.05,]
Cluster13_markers<-Cluster13_markers[order(Cluster13_markers$avg_log2FC,decreasing = T),]


Allcell_markers<-read.table("./epithelial_allcell_markers_bycluster.txt",sep="\t",header = T)
Featurename<-c("VAR01","VAR03","VAR05","VAR06","VAR07","Cluster13")
Featurename<-c("Bcell","Myeloid","Tcell","Plasma")
Features<-list()
# Features[[1]] <- Allcell_markers$gene[(Allcell_markers$cluster=="VAR01")&(Allcell_markers$avg_log2FC>=1)&(Allcell_markers$p_val_adj<0.05)]
# Features[[2]]<-Allcell_markers$gene[(Allcell_markers$cluster=="VAR03")&(Allcell_markers$avg_log2FC>=1)&(Allcell_markers$p_val_adj<0.05)]
# Features[[3]]<-Allcell_markers$gene[(Allcell_markers$cluster=="VAR05")&(Allcell_markers$avg_log2FC>=1)&(Allcell_markers$p_val_adj<0.05)]
# Features[[4]]<-Allcell_markers$gene[(Allcell_markers$cluster=="VAR06")&(Allcell_markers$avg_log2FC>=1)&(Allcell_markers$p_val_adj<0.05)]
# Features[[5]]<-Allcell_markers$gene[(Allcell_markers$cluster=="VAR07")&(Allcell_markers$avg_log2FC>=1)&(Allcell_markers$p_val_adj<0.05)]
# Features[[6]]<-Cluster13_markers$gene
Features[[1]]<-read.table("./Bcell_sig.csv",sep="\n",header = F)[[1]]
Features[[2]]<-read.table("./Myeloid_sig.csv",sep="\n",header = F)[[1]]
Features[[3]]<-read.table("./Tcell_sig.csv",sep="\n",header = F)[[1]]
Features[[4]]<-read.table("./Plasma_sig.csv",sep="\n",header = F)[[1]]
dge<-subset(dge,cells=colnames(dge)[!(dge@active.ident %in% c("VAR02","VAR08","UC03"))])
dge@active.ident<-factor()
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

pbmc.markers <- FindAllMarkers(dge,assay = "RNA",only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use = "MAST")
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
RPS.genes <- grep(pattern = "^RPS", x = rownames(dge@assays$RNA@data), value = TRUE)
RPL.genes <- grep(pattern = "^RPL", x = rownames(dge@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes,RPS.genes,RPL.genes,"MALAT1")
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
write.table(pbmc.markers,"DEG_byVAR.txt",col.names = T,row.names = T,quote = F,sep="\t")


top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.table(top10,"top10_markers_byVAR.txt",col.names = T,row.names = T,quote = F,sep="\t")


DotPlot(dge, features = rev(unique(top10$gene)), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave(file="Dotplot_10_byVAR.png",width = 15,height = 6)

dge@active.ident<-factor(dge@active.ident,levels = c("VAR01","VAR03","VAR05","VAR06","VAR07"))

dge@active.ident<-factor(dge@active.ident,levels = c("UC01","UC02","UC03","VAR01",
                                                     "VAR02","VAR03","VAR04","VAR05",
                                                     "VAR06","VAR07","VAR08","VAR09"))
dge@active.ident<-factor(dge@active.ident,levels = c("VAR09","VAR08","VAR07","VAR06",
                                                     "VAR05","VAR04","VAR03","VAR02",
                                                     "VAR01","UC03","UC02","UC01","Plasma"))
DotPlot(dge, features = c("Tcell1","Myeloid1","Plasma1","Bcell1"), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))+coord_flip()
ggsave("./Dotplot_Immune_with_plasma.pdf",height = 6,width = 8)


avg_expression<-AverageExpression(dge,return.seurat = T)
avg_expression<-as.matrix(avg_expression@assays$RNA@data)
avg_expression<-avg_expression[rownames(avg_expression) %in% Features[[4]],]
avg_expression<-as.data.frame(avg_expression)
avg_expression<-avg_expression[order(avg_expression$Plasma,decreasing = T),]
Top50_Plasma<-rownames(avg_expression)[1:50]
Top100_Plasma<-rownames(avg_expression)[1:100]
DotPlot(dge, features = Top50_Plasma, cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave("./Dotplot_Plasma.pdf",height = 6,width = 15)
ggsave("./Dotplot_Immune_flip.pdf",height = 6,width = 8)

DotPlot(dge, features = Top100_Plasma, cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave("./Dotplot_Plasma_top100.pdf",height = 6,width = 18)
###Project these back on the Epithelial cell object:
Epithelial<-readRDS("../Epithelial/dge_epithelial3_041023.rds")
Epithelial<-SetIdent(Epithelial,value = as.vector(Epithelial$Name2))
Epithelial@active.ident<-factor(Epithelial@active.ident,levels = c("UC01","UC02","UC03","VAR01","VAR02","VAR03","VAR04","VAR05","VAR06","VAR07","VAR08","VAR09","Cluster13"))
Featurename<-c("VAR01_Cluster13","VAR03_Cluster13","VAR05_Cluster13","VAR06_Cluster13","VAR07_Cluster13","Cluster13")
Features<-list()
Features[[1]] <- pbmc.markers$gene[(pbmc.markers$cluster=="VAR01")&(pbmc.markers$avg_log2FC>=0.5)&(pbmc.markers$p_val_adj<0.05)]
Features[[2]]<-pbmc.markers$gene[(pbmc.markers$cluster=="VAR03")&(pbmc.markers$avg_log2FC>=0.5)&(pbmc.markers$p_val_adj<0.05)]
Features[[3]]<-pbmc.markers$gene[(pbmc.markers$cluster=="VAR05")&(pbmc.markers$avg_log2FC>=0.5)&(pbmc.markers$p_val_adj<0.05)]
Features[[4]]<-pbmc.markers$gene[(pbmc.markers$cluster=="VAR06")&(pbmc.markers$avg_log2FC>=0.5)&(pbmc.markers$p_val_adj<0.05)]
#Features[[5]]<-pbmc.markers$gene[(pbmc.markers$cluster=="VAR07")&(pbmc.markers$avg_log2FC>=0.5)&(pbmc.markers$p_val_adj<0.05)]
Features[[5]]<-top20$gene[top20$cluster=="VAR07"]
Features[[6]]<-Cluster13_markers$gene

for (i in 1:length(Features)) {
  gene.set <- Features[[i]][1:length(Features[[i]])]
  Epithelial <- AddModuleScore(object = Epithelial, features = list(gene.set), ctrl = 2, name = Featurename[i])
  V1<-as.vector(FetchData(object = Epithelial,vars = paste0(Featurename[i],"1")))
  FeaturePlot(Epithelial,features=paste0(c(Featurename[i]),"1"))+scale_color_gradientn( colours = c("blue","green","yellow","red"),  limits = c(0, max(unlist(V1))))
  ggsave(file=paste0("Featureplot_",Featurename[i],".pdf"),width = 8,height = 8)
  BOX_df<-NULL
  BOX_df$id<-Epithelial@active.ident
  BOX_df$value<-as.vector(V1)
  BOX_df<-data.frame(BOX_df)
  colnames(BOX_df)<-c("id","value")
  ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot() + #geom_violin() +#geom_jitter(shape=16,position = position_jitter(0.1))+
    stat_summary(fun=mean,geom="point",size=3,colour="blue",shape=95)+ theme(legend.position="none",text = element_text(size=6))+
    stat_compare_means(method = "anova",label.x = 3,label.y = max(unlist(V1))+0.05)+
    ggtitle(Featurename[i])+rotate_x_text(angle = 45)
  ggsave(file=paste0("Vlnplot_",Featurename[i],".pdf"),width = 10,height = 5)
  
}
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




####ENhances volcano
setwd("~/Desktop/Bladder_2023/")
res<-read.table("./Volcano_dge_variantvsUC.pdf_Variant_Pure UC_Volcano.txt",sep="\t",header = T)
res<-res[res$cluster=="Variant",]
rownames(res)<-res$gene
source("~/Desktop/GO_loose_analysis_geneset.R")
res<-res[order(res$avg_log2FC,decreasing = T),]
Markers<-res
name_target<-"Variant"
threshold=1
threshold_2=1
p_threshold=0.05
GO_loose_analysis(Markers,name_target,threshold,threshold_2,p_threshold)

colnames(res)[2]="log2FoldChange"
colnames(res)[5]="pvalue"
res<-res[order(res$cluster,decreasing = T),]
#res<-res[1:(nrow(res)/2),]
#res$pvalue[res$pvalue==0]=1.689000e-301
library(EnhancedVolcano)
#res$pvalue[res$pvalue==0]=5.944247e-303
#res$pvalue[res$pvalue==0]=6.813668e-296
options(ggrepel.max.overlaps = Inf)
Genes<-res[order(abs(res$log2FoldChange),decreasing = T),]
Genes<-Genes$gene[1:20]
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = Genes,
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-10,
                FCcutoff = 2.0,
                pointSize = 0.5,
                labSize = 3.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                maxoverlapsConnectors = Inf,
                axisLabSize = 10)


EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c(Genes),
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-10,
                FCcutoff = 2.0,
                pointSize = 0.5,
                labSize = 2,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.2,
                colConnectors = 'black',
                maxoverlapsConnectors = Inf,
                axisLabSize = 10) +NoLegend()


###Correlation
dge<-readRDS("/Volumes/Backup_Plus/Bladder_2023/Epithelial.rds")
HLA_class1<-c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G")
HLA_class2<-c("HLA-DOA","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","HLA-DMA","HLA-DMB","HLA-DOB","HLA-DQA2","HLA-DQB2","HLA-F-AS1","HLA-DQB1-AS1","HLA-AS1")

dge <- AddModuleScore(object = dge, features = list(HLA_class1), ctrl = 2, name = "HLA_classI")
dge <- AddModuleScore(object = dge, features = list(HLA_class2), ctrl = 2, name = "HLA_classII")

Gene1="TM4SF1"
Gene2="CD274" ###change it to HLA_classII1
Expression_Gene1<-FetchData(dge,vars = Gene1)
Expression_Gene2<-FetchData(dge,vars = Gene2)
Meta<-dge$Name
#cor.test(as.numeric(Expression_Gene1[,1]),as.numeric(Expression_Gene2[,1]))
DF<-cbind(Expression_Gene1,Expression_Gene2,Meta)
colnames(DF)<-c("Gene1","Gene2","Group")
DF$Gene1<-as.numeric(DF$Gene1)
DF$Gene2<-as.numeric(DF$Gene2)
ggplot(DF, aes(y = Gene2, x = Gene1)) +
  geom_point(aes(colour = factor(Group))) +
  facet_wrap("Group") +
  geom_smooth(method=lm) + 
  xlab(Gene1) + ylab(Gene2) +NoLegend()
ggsave(paste0("Correlation_",Gene2,"_",Gene1,".pdf"),width = 8,height = 6)
###

###Compute correlation with TM4SF1
library("Hmisc")
library(corrplot)
dge<-readRDS("./dge_epithelial_032823.rds")
dge<-SetIdent(dge,value = as.vector(dge$Broad_type))
ID="Pure UC"
dge<-subset(dge,idents=c(ID))
dge <- NormalizeData(object = dge, normalization.method = "LogNormalize", scale.factor = 10000)
dge <- ScaleData(object=dge,features=rownames(dge))

matrix<-dge@assays$RNA@data
matrix_mod<-as.matrix(matrix)

gene<-as.numeric(matrix_mod["FTH1",])
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
write.table(DF,paste0("Correlation_FTH1.txt"),sep="\t",col.names = T,row.names = T,quote = F)


Tumor_DF<-read.table("./Correlation_TM4SF1_Variant.txt",sep="\t",header = T)
Tumor_DF<-Tumor_DF[Tumor_DF$P_value<0.05,]
Tumor_DF<-Tumor_DF[abs(Tumor_DF$Correlation)>=0.2,]
Tumor_DF<-Tumor_DF[!is.na(Tumor_DF$P_value),]
Tumor_DF<-Tumor_DF[order(Tumor_DF$Correlation,decreasing = T),]
Tumor_DF<-Tumor_DF[-c(1),]

NonTumor_DF<-read.table("./Correlation_TM4SF1_Pure UC.txt",sep="\t",header = T)
# NonTumor_DF<-NonTumor_DF[NonTumor_DF$P_value<0.05,]
# NonTumor_DF<-NonTumor_DF[NonTumor_DF$Correlation>=0.2,]
# NonTumor_DF<-NonTumor_DF[!is.na(NonTumor_DF$P_value),]
# NonTumor_DF<-NonTumor_DF[order(NonTumor_DF$Correlation,decreasing = T),]
# NonTumor_DF<-NonTumor_DF[-c(1),]

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
write.table(DF,"./TM4SF1_Correlation.txt",sep = "\t",col.names = T,row.names = F,quote = F)
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

####Luminal_Basal
Luminal<-c("UPK2","UPK1A","UPK1B")
Basal<-c("KRT5","KRT6A","KRT15")
Luminal_Guo <- c("CYP2J2","ERBB2","ERBB3","FGFR3","FOXA1","GATA3","GPX2","KRT18","KRT19","KRT20","KRT7","KRT8","PPARG","XBP1")
Basal_Guo <- c("CD44","CDH3","KRT1","KRT14","KRT16","KRT5","KRT6A","KRT6B","KRT6C")

dge <- AddModuleScore(object = dge, features = list(Luminal), ctrl = 2, name = "Luminal")
dge <- AddModuleScore(object = dge, features = list(Basal), ctrl = 2, name = "Basal")
dge <- AddModuleScore(object = dge, features = list(Luminal_Guo), ctrl = 2, name = "Luminal_Guo")
dge <- AddModuleScore(object = dge, features = list(Basal_Guo), ctrl = 2, name = "Basal_Guo")


dge$Luminal_Basal="Others"
dge$Luminal_Basal[(dge$Luminal1>0)&(dge$Basal1<=0)]="Luminal-like"
dge$Luminal_Basal[(dge$Luminal1<=0)&(dge$Basal1>0)]="Basal-like"

dge$all = "all"
Gene1="TM4SF1"
Gene2="CLDN4"
Expression_Gene1<-FetchData(dge,vars = Gene1)
Expression_Gene2<-FetchData(dge,vars = Gene2)
Meta<-dge$Name2
#cor.test(as.numeric(Expression_Gene1[,1]),as.numeric(Expression_Gene2[,1]))
DF<-cbind(Expression_Gene1,Expression_Gene2,Meta)
colnames(DF)<-c("Gene1","Gene2","Group")
DF$Gene1<-as.numeric(DF$Gene1)
DF$Gene2<-as.numeric(DF$Gene2)
ggplot(DF, aes(y = Gene2, x = Gene1)) +
  geom_point(aes(colour = factor(Group)),size=0.5) +
  facet_wrap("Group") +
  geom_smooth(method=lm) +
  xlab(Gene1) + ylab(Gene2) +theme(panel.background = element_rect(fill = "white")) + NoLegend()
ggsave(paste0("Correlation_",Gene1,"_",Gene2,"_Name2.pdf"),width = 10,height = 7)


####TM4SF1
DF_variant<-read.table("./Correlation_FTH1.txt",sep="\t",header = T,row.names = 1)
DF_variant<-DF_variant[!is.na(DF_variant$Correlation),]
DF_variant$FDR<-p.adjust(DF_variant$P_value,n = nrow(DF_variant), method = "fdr")
# DF_variant$FDR[DF_variant$FDR==0]=1.762287e-15
# DF_variant$log10FDR<-(-1)*(log10(DF_variant$FDR))
DF_variant<-DF_variant[order(DF_variant$Correlation,decreasing = T),]
DF_variant<-DF_variant[-c(1),]
DF_variant$label=0
DF_variant$label[1:10]=1
DF_variant$label[(nrow(DF_variant)-9):nrow(DF_variant)]=1
ggplot(DF_variant, aes(x = Correlation, y = log10FDR)) +
  geom_point(aes(color = ifelse(label==1, "highlight", "normal")),size=0.5) +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) + 
  geom_label_repel(aes(label = ifelse(label==1, Gene, "")), 
                  # fill = "white",label.padding = 0,box.padding = 0,
                   size = 2,
                   fontface = "bold",
                   segment.color = "grey50",
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
                   segment.size = 0.2,hjust = ifelse(DF_variant$label==1, 0, 1),
                   max.overlaps = 50) + NoLegend()

Variant<-read.table("./Correlation_TM4SF1_Variant.txt",sep="\t",header = T,row.names = 1)
UC<-read.table("./Correlation_TM4SF1_Pure UC.txt",sep="\t",header = T,row.names = 1)
UC<-UC[!is.na(UC$P_value),]
Variant<-Variant[!is.na(Variant$P_value),]
DF<-NULL
Genelist<-intersect(Variant$Gene,UC$Gene)
DF$Gene<-Genelist
DF$Variant<-0
DF$UC<-0
for (i in 1:length(Genelist)){
  DF$Variant[i]=Variant$Correlation[Variant$Gene==Genelist[i]]
  DF$UC[i]=UC$Correlation[UC$Gene==Genelist[i]]
}
DF<-as.data.frame(DF)
Variant<-Variant[order(Variant$Correlation,decreasing = T),]
Variant<-Variant[-c(1),]
Gene_highlight_variant<-Variant$Gene[1:10]

UC<-UC[order(UC$Correlation,decreasing = T),]
UC<-UC[-c(1),]
Gene_highlight_UC<-UC$Gene[1:10]

DF$label=0
DF$label[DF$Gene %in% Gene_highlight_variant]=1
DF$label[DF$Gene %in% Gene_highlight_UC]=1
ggplot(DF, aes(x = Variant, y = UC)) +
  geom_point(aes(color = ifelse(label==1, "highlight", "normal")),size=0.5) +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) + 
  geom_label_repel(aes(label = ifelse(label==1, Gene, "")), 
                   # fill = "white",label.padding = 0,box.padding = 0,
                   size = 2,
                   fontface = "bold",
                   segment.color = "grey50",
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
                   segment.size = 0.2,#hjust = ifelse(DF$label==1, 0, 1),
                   max.overlaps = 50) + NoLegend()


####scatterplot of correlation vs expression
dge<-SetIdent(dge,value = as.vector(dge$ID))
Avg<-AverageExpression(dge)
Avg<-as.data.frame(Avg)
DF<-merge(Avg,DF,by.x="row.names",by.y="Gene")
colnames(DF)<-c("Gene","Expression","Correlation","P_value")
DF<-DF[!is.na(DF$P_value),]
DF<-DF[order(DF$Correlation,decreasing = T),]
DF<-DF[-c(1),]
mito.genes <- grep(pattern = "^MT-", x = rownames(dge@assays$RNA@data), value = TRUE)
DF<-DF[!(DF$Gene %in% c(mito.genes,"MALAT1","RNA28S5")),]
DF$label=0
DF$label[1:10]=1
DF$label[((nrow(DF)-9):nrow(DF))]=1
DF$label[DF$Gene %in% c("KLK4","MIF")]=1
ggplot(DF, aes(x = Correlation, y = Expression)) +
  geom_point(aes(color = ifelse(label==1, "highlight", "normal")),size=0.5) +
  scale_color_manual(values = c("highlight" = "red", "normal" = "black")) + 
  geom_label_repel(aes(label = ifelse(label==1, Gene, "")), 
                   # fill = "white",label.padding = 0,box.padding = 0,
                   size = 2,
                   fontface = "bold",
                   segment.color = "grey50",
                   arrow = arrow(length = unit(0.01, "npc"), type = "closed"),
                   segment.size = 0.2,#hjust = ifelse(DF$label==1, 0, 1),
                   max.overlaps = 50) + 
  theme(axis.text=element_text(size=8),axis.title=element_text(size=8)) + 
  NoLegend()

geneset<-DF$Gene[DF$Correlation>0.25]

####CIBERSORT
DF<-DF[DF$cis %in% c("0","1"),]

ggplot(data=DF, aes(x=Cluster13, y=cis)) +
  geom_point(aes(color=cis), size=5) + 
  geom_smooth(method = "glm", se = T, method.args = list(family = "binomial"), linetype = "dashed") + 
  xlab("Cluster13_Percentage") +
  ylab("CIS_Status")
ggsave("CIS_Cluster13_GLM.pdf")


ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot()

DF$cis<-as.character(DF$cis)
ggplot(DF, aes(cis,Cluster13,fill=cis))  + geom_boxplot() + facet_wrap("variant_category", ncol = 4) + geom_jitter(shape=16,position = position_jitter(0.1)) +
  stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95) +
  geom_signif(comparisons = list(c("0","1")),map_signif_level=TRUE,y_position = max(DF$Cluster13), test = "wilcox.test",alternative = "two.sided",paired = T) +
  #geom_text(data=DF_feature, hjust=0, size=3, aes(x=1,y=round(max(log2count, na.rm=TRUE)+0.5), label=paste0("Wilcox test p=",normData$pvalue[normData$GeneSymbol==Featurename]))) + 
  ggtitle("CIS_Cluster13") + NoLegend()+ theme_bw()


