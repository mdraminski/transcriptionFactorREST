
## ---------------------------

## Purpose of script: Comparison of IDH-MUT and IDH-WT in TCGA gliomas and U87 cell model - comparison on the level of REACTOME Pathways
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

library(DESeq2)

load("/media/bartek/OS/TCGA_ATRX_for_MCFS/TCGA_LGG_GBM_DESeq2_own_normalization_IDH_MUT_vs_WT_DESeq2_analysis_DEC_2021_plus_G4_WT_fr3.RData")

##Read in libraries
library("clusterProfiler")
library("ggplot2")
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)

w1<-which(fr1_lim$padj<0.01 & fr1_lim$log2FoldChange>0)
w2<-which(fr2_lim$padj<0.01 & fr2_lim$log2FoldChange>0)
w3<-which(fr3_lim$padj<0.01 & fr3_lim$log2FoldChange>0)
m1<-match(rownames(fr1_lim)[w1],gen$ensembl_gene_id)
m2<-match(rownames(fr2_lim)[w2],gen$ensembl_gene_id)
m3<-match(rownames(fr3_lim)[w3],gen$ensembl_gene_id)

# Create a list with genes from each sample
genes = list(na.omit(gen[m1,9]),na.omit(gen[m2,9]))
names(genes)<-c("IDH_WT","U87_MUT_WT")

compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "fdr")
dotplot(compKEGG, showCategory = 20, title = "REACTOME Pathway Enrichment Analysis")




w1<-which(fr1_lim$padj<0.01 & fr1_lim$log2FoldChange< 0)
w2<-which(fr2_lim$padj<0.01 & fr2_lim$log2FoldChange< 0)
w3<-which(fr3_lim$padj<0.01 & fr3_lim$log2FoldChange< 0)
m1<-match(rownames(fr1_lim)[w1],gen$ensembl_gene_id)
m2<-match(rownames(fr2_lim)[w2],gen$ensembl_gene_id)
m3<-match(rownames(fr3_lim)[w3],gen$ensembl_gene_id)

# Create a list with genes from each sample
genes = list(as.character(na.omit(gen[m1,9])),as.character(na.omit(gen[m2,9])))
names(genes)<-c("IDH_WT","U87_MUT_WT")

compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichPathway",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "fdr")                        
dotplot(compKEGG, showCategory = 20, title = "REACTOME Pathway Enrichment Analysis")


