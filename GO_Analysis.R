
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
hallmark<-read.gmt("~/Desktop/Geneset/hallmark.gmt")
c2_cp<-read.gmt("~/Desktop/Geneset/c2.cp.v7.4.symbols.gmt")
c2_kegg<-read.gmt("~/Desktop/Geneset/c2.KEGG.v7.4.symbols.gmt")
c2_cgp<-read.gmt("~/Desktop/Geneset/c2.cgp.v7.4.symbols.gmt")
c4<-read.gmt("~/Desktop/Geneset/c4.all.v7.4.symbols.gmt")
c6<-read.gmt("~/Desktop/Geneset/c6.all.v7.4.symbols.gmt")

Name="Outlier"
genes<-read.table("./Outlier_DEG.txt",sep="\n",header = F)
sig_cluster<-genes
gene_symbols<-as.vector(t(sig_cluster))
ego_symbols <- enrichGO(gene_symbols, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont = "ALL", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)
dotplot(ego_symbols,showCategory=20,font.size=25)+ggtitle("C5 Ontology")
ggsave(paste0(Name,"_Enrichment_C5.pdf"),height = 20,width = 10)

egmt <- enricher(gene_symbols, TERM2GENE=hallmark,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dotplot(egmt,showCategory=20,font.size=20)+ggtitle("Hallmark")
ggsave(paste0(Name,"_Enrichment_Hallmark.pdf"),height = 10,width = 8)

egmt <- enricher(gene_symbols, TERM2GENE=c2_cp,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dotplot(egmt,showCategory=20,font.size=25)+ggtitle("C2CP")
ggsave(paste0(Name,"_Enrichment_C2CP.pdf"),height = 20,width = 12)

egmt <- enricher(gene_symbols, TERM2GENE=c2_kegg,pAdjustMethod = "BH",qvalueCutoff = 0.05)
dotplot(egmt,showCategory=20,font.size=25)+ggtitle("C2KEGG")
ggsave(paste0(Name,"_Enrichment_C2CPKEGG.pdf"),height = 15,width = 12)

egmt <- enricher(gene_symbols, TERM2GENE=c2_cgp)
dotplot(egmt,showCategory=20,font.size=25)+ggtitle("C2CGP")
ggsave(paste0(Name,"_Enrichment_C2CGP.pdf"),height = 15,width = 12)
