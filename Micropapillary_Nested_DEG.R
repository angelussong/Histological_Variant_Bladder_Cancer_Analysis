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

#####Investigate Common Demoninators of Nested or Micropapillary 
dge<-readRDS("./dge_epithelial2final.rds")
mkdir("./Common_Demoninator")
setwd("./Common_Demoninator/")
dge<-SetIdent(dge,value = as.vector(dge$histology2))
PureUC_ID<-c("Pure Urothelial (12923)","Pure Urothelial (SG)","Pure Urothelial (HG)")
Target_ID<-c("Micropapillary+pleomorphic (JM)")
dge_subset<-subset(dge,idents=c(PureUC_ID,Target_ID))
dge_subset$Target="PureUC"
dge_subset$Target[dge_subset@active.ident %in% Target_ID]=Target_ID
dge_subset<-SetIdent(dge_subset,value = as.vector(dge_subset$Target))
pbmc.markers <- FindAllMarkers(dge_subset,assay = "RNA",only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
mito.genes <- grep(pattern = "^MT-", x = rownames(dge_subset@assays$RNA@data), value = TRUE)
remove_gene<-c(mito.genes)
pbmc.markers<-pbmc.markers[!(pbmc.markers$gene %in% remove_gene),]
write.table(pbmc.markers,paste0(Target_ID,"_Markers.txt"),col.names = T,row.names = T,quote = F,sep="\t")

top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
markers.to.plot <- unique(top10$gene)
DotPlot(dge_subset, features = unique(markers.to.plot), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave(file=paste0(Target_ID,"_Markers.pdf"),width = 30,height = 10,units = "cm",limitsize = FALSE)

# Nested_FG_markers<-pbmc.markers[pbmc.markers$cluster==Target_ID,]
# Nested_FG_markers<-Nested_FG_markers[Nested_FG_markers$p_val_adj<0.05,]
# Nested_FG_markers<-Nested_FG_markers[Nested_FG_markers$avg_log2FC>0.5,]
# 
# 
# Nested_12041_markers<-pbmc.markers[pbmc.markers$cluster==Target_ID,]
# Nested_12041_markers<-Nested_12041_markers[Nested_12041_markers$p_val_adj<0.05,]
# Nested_12041_markers<-Nested_12041_markers[Nested_12041_markers$avg_log2FC>0.5,]

Micro_PS_markers<-pbmc.markers[pbmc.markers$cluster==Target_ID,]
Micro_PS_markers<-Micro_PS_markers[Micro_PS_markers$p_val_adj<0.05,]
Micro_PS_markers<-Micro_PS_markers[Micro_PS_markers$avg_log2FC>0.5,]

Micro_21222_markers<-pbmc.markers[pbmc.markers$cluster==Target_ID,]
Micro_21222_markers<-Micro_21222_markers[Micro_21222_markers$p_val_adj<0.05,]
Micro_21222_markers<-Micro_21222_markers[Micro_21222_markers$avg_log2FC>0.5,]

Micro_JM_markers<-pbmc.markers[pbmc.markers$cluster==Target_ID,]
Micro_JM_markers<-Micro_JM_markers[Micro_JM_markers$p_val_adj<0.05,]
Micro_JM_markers<-Micro_JM_markers[Micro_JM_markers$avg_log2FC>0.5,]

# Common_Nested<-intersect(Nested_12041_markers$gene,Nested_FG_markers$gene)
# dge_Target<-subset(dge,idents=c(PureUC_ID,"Nested (12041)","Nested (FG)"))
# dge_Target$Target=as.vector(dge_Target@active.ident)
# dge_Target$Target[dge_Target@active.ident %in% PureUC_ID]="PureUC"
# dge_Target<-SetIdent(dge_Target,value = as.vector(dge_Target$Target))
# DotPlot(dge_Target, features = unique(Common_Nested), cols = "RdBu")+
#   theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
# ggsave(file=paste0("Common_Nested","_Markers.pdf"),width = 40,height = 10,units = "cm",limitsize = FALSE)
# 
# write.table(Common_Nested,"Common_Nested_Markers.txt",sep="\n",col.names = F,row.names = F,quote = F)

Common_Micro<-intersect(Micro_PS_markers$gene,Micro_21222_markers$gene)
Common_Micro_three<-intersect(Common_Micro,Micro_JM_markers$gene)
dge_Target<-subset(dge,idents=c(PureUC_ID,"Micropapillary (21222)","Micropapillary (PS)","Micropapillary+pleomorphic (JM)"))
dge_Target$Target=as.vector(dge_Target@active.ident)
dge_Target$Target[dge_Target@active.ident %in% PureUC_ID]="PureUC"
dge_Target<-SetIdent(dge_Target,value = as.vector(dge_Target$Target))
DotPlot(dge_Target, features = unique(Common_Micro_three), cols = "RdBu")+
  theme(axis.text.x = element_text( angle = 90, vjust = 0.5, hjust = 1))
ggsave(file=paste0("Common_Micro","_Markers.pdf"),width = 70,height = 10,units = "cm",limitsize = FALSE)
ggsave(file=paste0("Common_Micro","_Markers_three.pdf"),width = 30,height = 10,units = "cm",limitsize = FALSE)

write.table(Common_Micro_three,"Common_Micro_Markers_three.txt",sep="\n",col.names = F,row.names = F,quote = F)

