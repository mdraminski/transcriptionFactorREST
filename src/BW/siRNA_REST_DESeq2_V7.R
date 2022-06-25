
## ---------------------------

## Purpose of script: Analysis of siREST effect on U87 ID-MUT and IDH-WT glioma cell lines
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_siRNA")
counts<-read.table("feature_counts_siRNA_REST_2020", sep="\t",fill=TRUE,header=TRUE,skip=1,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
data<-counts[grep("ENST",counts$Geneid,invert = TRUE),c(7:ncol(counts))]
rownames(data)<-gsub("gene:","",counts[grep("ENST",counts$Geneid,invert = TRUE),1])
a<-gsub(".bam","",colnames(data))
colnames(data)<-a

WT_siCTRL<-intersect(grep("F11",colnames(data)),grep("_V11",colnames(data)))
WT_siREST<-intersect(grep("F11",colnames(data)),grep("_V21",colnames(data)))
MUT_siCTRL<-intersect(grep("F31",colnames(data)),grep("_V11",colnames(data)))
MUT_siREST<-intersect(grep("F31",colnames(data)),grep("_V21",colnames(data)))

mycounts<-cbind(data[,WT_siCTRL],data[,WT_siREST],data[,MUT_siCTRL],data[,MUT_siREST])

###Annotating 
library(biomaRt)
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='jul2018.archive.ensembl.org')
gen<- getBM(values = rownames(mycounts), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description","entrezgene"))


library("DESeq2")
condition <- factor(c(rep("WT_siCTRL",3),rep("WT_siREST",3),rep("MUT_siCTRL",3),rep("MUT_siREST",3)))
Patient.ID <- factor(c(1,2,3,1,2,3,1,2,3,1,2,3))


#coldata <- data.frame(row.names=colnames(mycounts), condition)
coldata <- data.frame(row.names=colnames(mycounts), condition, Patient.ID)

#create a DESeq object (New:multifactoral design including the patients ID to account for differences between the samples
#Fit an individual base line for each patient, so that patient-to-patient differences in expression before treatment are absorbed and you only look at the differences the condition)              
dds1 <- DESeqDataSetFromMatrix(countData=mycounts, colData=coldata, design=~Patient.ID + condition)
dds1
dds1 <- DESeq(dds1)
fr1<-results(dds1,contrast = c("condition", "WT_siREST","WT_siCTRL"))
fr2<-results(dds1,contrast = c("condition", "MUT_siREST","MUT_siCTRL"))
fr3<-results(dds1,contrast = c("condition","MUT_siREST","WT_siREST"))
fr4<-results(dds1,contrast = c("condition","MUT_siCTRL","WT_siCTRL"))
###hist(fr4$padj) and hist(fr3$padj) are almost identical - showing differences U87-MUT vs U87-WT


setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_siRNA/2ndary")
m<-match(fr1@rownames,gen[,1])
res<-cbind(fr1@rownames,gen[m,],fr1$log2FoldChange,fr1$pvalue,fr1$padj,
           fr2$log2FoldChange,fr2$pvalue,fr2$padj,
           fr3$log2FoldChange,fr3$pvalue,fr3$padj,
           fr4$log2FoldChange,fr4$pvalue,fr4$padj)
colnames(res)<-c("ens",colnames(gen),"log2FoldChange_WT_siREST_vs_WT_siCTRL","pvalue_WT_siREST_vs_WT_siCTRL","padj_WT_siREST_vs_WT_siCTRL",
                 "log2FoldChange_MUT_siREST_vs_MUT_siCTRL","pvalue_MUT_siREST_vs_MUT_siCTRL","padj_MUT_siREST_vs_MUT_siCTRL",
                 "log2FoldChange_MUT_siREST_vs_WT_siREST","pvalue_MUT_siREST_vs_WT_siREST","padj_MUT_siREST_vs_WT_siREST",
                 "log2FoldChange_MUT_siCTRL_vs_WT_siCTRL","pvalue_MUT_siCTRL_vs_WT_siCTRL","padj_MUT_siCTRL_vs_WT_siCTRL")




###get matrix of normalized counts
df <- as.data.frame(counts(dds1, normalized=TRUE))
m<-match(rownames(df),gen[,1])


library(RColorBrewer)
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)


####Improved DEC 20202

#####Volcano plot: WT siREST/siCTRL
res2 <- results(dds1,contrast = c("condition","WT_siREST", "WT_siCTRL"))
m<-match(rownames(res2),gen[,1])
res3<-res2[which(!is.na(m)),]
rownames(res3)<-gen[na.omit(m),2]

threshold_OE <- rep(NA,nrow(res3))
w1<-which(res3$padj<0.05 & res3$log2FoldChange>0)
w2<-which(res3$padj<0.05 & res3$log2FoldChange<0)
threshold_OE[w1]=Colors[100]
threshold_OE[w2]=Colors[1]
threshold_OE[-c(w1,w2)]="lightgrey"


res3$threshold <- threshold_OE
results_ordered <- res3[order(res3$padj), ]

results_ordered$genelabels <- ""
results_ordered$genelabels[1:nrow(subset(results_ordered, c(w1,w2)))] <- rownames(results_ordered)[1:nrow(subset(results_ordered, c(w1,w2)))]
#Volcano plot
library(ggrepel)
title = ""
xtitle = "log2 Fold change siREST/siCTRL"
#x_limits <- c(3, NA)

mycolors <- c(Colors[100],Colors[1], "#797b77")
names(mycolors) <- c(Colors[100],Colors[1], "#797b77")


pdf(file="/home/bartek/Pulpit/podsumowanie_RES_V7_ALL_FIGURES/volcano_siREST_siCTRL_WT_V7.pdf", height=10, width=6)
ggplot(as.data.frame(results_ordered)) +
  theme_bw() +
  scale_color_manual(values = mycolors) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(size=7,segment.size = 0.5,nudge_y = 0.1, force = 10, segment.alpha = 0.5, direction = "both",hjust=0,box.padding = 0.5, aes(x = log2FoldChange, y = -log10(padj),label = ifelse((results_ordered$genelabels != "" & results_ordered$log2FoldChange < -1.5 & results_ordered$padj < 0.05), genelabels,""))) +
  geom_text_repel(size=7,segment.size = 0.5,nudge_y = 0.1, force = 10, segment.alpha = 0.5, direction = "both",hjust=1,box.padding = 0.5, aes(x = log2FoldChange, y = -log10(padj),label = ifelse((results_ordered$genelabels != "" & results_ordered$log2FoldChange > 2 & results_ordered$padj < 0.05), genelabels,""))) +
  ggtitle(title) +
  xlab(xtitle) +
  ylab("Statistical significance (-log10 q-value)") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(2))) + xlim(-2.5,8) +
  theme(axis.text.x= element_text(size=14)) +
  theme(axis.text.y= element_text(size=14)) +
  theme(
    panel.border = element_blank()
  )
dev.off()




#####Volcano plot: MUT siREST/siCTRL
res2 <- results(dds1,contrast = c("condition","MUT_siREST", "MUT_siCTRL"))
m<-match(rownames(res2),gen[,1])
res3<-res2[which(!is.na(m)),]
rownames(res3)<-gen[na.omit(m),2]

#Volcano plot - Adria's code

threshold_OE <- rep(NA,nrow(res3))
w1<-which(res3$padj<0.05 & res3$log2FoldChange>0)
w2<-which(res3$padj<0.05 & res3$log2FoldChange<0)
threshold_OE[w1]=Colors[100]
threshold_OE[w2]=Colors[1]
threshold_OE[-c(w1,w2)]="lightgrey"


res3$threshold <- threshold_OE
results_ordered <- res3[order(res3$padj), ]

results_ordered$genelabels <- ""
results_ordered$genelabels[1:nrow(subset(results_ordered, c(w1,w2)))] <- rownames(results_ordered)[1:nrow(subset(results_ordered, c(w1,w2)))]
#Volcano plot
library(ggrepel)
title = ""
xtitle = "log2 Fold change siREST/siCTRL"
#x_limits <- c(3, NA)

mycolors <- c(Colors[100],Colors[1], "#797b77")
names(mycolors) <- c(Colors[100],Colors[1], "#797b77")


pdf(file="/home/bartek/Pulpit/podsumowanie_RES_V7_ALL_FIGURES/volcano_siREST_siCTRL_MUT_V7.pdf", height=10, width=6)
ggplot(as.data.frame(results_ordered)) +
  theme_bw() +
  scale_color_manual(values = mycolors) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold)) +
  geom_text_repel(size=7,segment.size = 0.5,nudge_y = 0.1, force = 10, segment.alpha = 0.5, direction = "both",hjust=0,box.padding = 0.5, aes(x = log2FoldChange, y = -log10(padj),label = ifelse((results_ordered$genelabels != "" & results_ordered$log2FoldChange < -1.5 & results_ordered$padj < 0.05), genelabels,""))) +
  geom_text_repel(size=7,segment.size = 0.5,nudge_y = 0.1, force = 10, segment.alpha = 0.5, direction = "both",hjust=1,box.padding = 0.5, aes(x = log2FoldChange, y = -log10(padj),label = ifelse((results_ordered$genelabels != "" & results_ordered$log2FoldChange > 2 & results_ordered$padj < 0.05), genelabels,""))) +
  ggtitle(title) +
  xlab(xtitle) +
  ylab("Statistical significance (-log10 q-value)") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(2), hjust = 0.5),
        axis.title = element_text(size = rel(2))) + xlim(-2.5,8) +
  theme(axis.text.x= element_text(size=14)) +
  theme(axis.text.y= element_text(size=14)) +
  theme(
    panel.border = element_blank()
  )
dev.off()


library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


w1<-which(fr1$padj<0.05 & fr2$padj<0.05 & fr1$log2FoldChange>0 & fr2$log2FoldChange>0)
w2<-which(fr1$padj<0.05 & fr2$padj<0.05 & fr1$log2FoldChange<0 & fr2$log2FoldChange<0)

m<-match(fr1@rownames[w1],gen[,1])

kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Up-regulated siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))

m<-match(fr1@rownames[w2],gen[,1])
kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Down-regulated siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))







#####additional for Gosia


w1<-which(fr1$padj<0.05 & fr1$log2FoldChange>0 )
w2<-which(fr1$padj<0.05 & fr1$log2FoldChange<0 )

m<-match(fr1@rownames[w1],gen[,1])

kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Up-regulated WT siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))

m<-match(fr1@rownames[w2],gen[,1])

kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Down-regulated WT siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))




w1<-which(fr2$padj<0.05 & fr2$log2FoldChange>0 )
w2<-which(fr2$padj<0.05 & fr2$log2FoldChange<0 )

m<-match(fr2@rownames[w1],gen[,1])

kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Up-regulated MUT siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))

m<-match(fr2@rownames[w2],gen[,1])

kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb = "org.Hs.eg.db",
                         pvalueCutoff = 0.05,ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO Down-regulated MUT siREST vs siCTRL") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))

