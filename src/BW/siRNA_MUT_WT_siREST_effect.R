
## ---------------------------

## Purpose of script: Comparison of IDH-MUT and IDH-WT in not-treaded, siCTRl- or siREST-treated IDH-MUT vs IDH-WT
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

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


####WT_siREST_vs_WT_siCTRL

w1<-which(fr3$padj<.05 & (fr3$log2FoldChange-fr4$log2FoldChange)< -0.25)
w2<-which(fr3$padj<.05 & (fr3$log2FoldChange-fr4$log2FoldChange)>0.25)
w3<-which(fr3$padj<.05 & (fr3$log2FoldChange-fr4$log2FoldChange)<0.25 & (fr3$log2FoldChange-fr4$log2FoldChange)> -0.25)
            

# 3D Exploded Pie Chart
library(plotrix)
slices <- c(length(w3), length(w1),length(w2))
lbls <- c("", "", "")
pie3D(slices,labels=lbls,explode=0.25,col=c("lightgrey","red","green"),
      main="fold changes MUT vs WT",mar=c(1,1,1,1))

setwd("/home/bartek/Pulpit/podsumowanie_RES_V7_ALL_FIGURES")



m<-match(fr3@rownames[w1],gen[,1])
kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb="org.Hs.eg.db",
               pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO BP Lower in MUT when REST is silenced") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))


GENELIST=fr3$log2FoldChange-fr4$log2FoldChange
names(GENELIST)<-res$entrezgene

w<-which(fr1$padj<0.05)
edox <- setReadable(kk_nfatc1_up, 'org.Hs.eg.db', 'ENTREZID')
###Wersja "gałązka"
cnetplot(edox, foldChange=GENELIST)
cnetplot(edox, foldChange=GENELIST) +
  scale_color_gradient(low = "#ff0000", high = "#ffa500",name="custom")


par(mar=c(6,6,1,1))
plot(fr4$log2FoldChange[w3],fr3$log2FoldChange[w3], col="lightgrey",pch=19, ylab="log2Change_MUT_WT_siREST",xlab="log2Change_MUT_WT_siCTRL",
     main="",cex.axis=2,cex.lab=2,bty="n")
points(fr4$log2FoldChange[w2],fr3$log2FoldChange[w2], col="green",pch=19)
points(fr4$log2FoldChange[w1],fr3$log2FoldChange[w1], col="red",pch=19)
legend("bottomright",legend = c("siREST independent","iDEGs","dDEGs"),col=c("lightgrey","green","red"),pch=19,bty='n',cex=1.5)




m<-match(fr3@rownames[w2],gen[,1])
kk_nfatc1_up <- enrichGO(gene         = gen$entrezgene[m],OrgDb="org.Hs.eg.db",
                         pvalueCutoff = 0.05,pAdjustMethod="fdr",ont="BP")
barplot(kk_nfatc1_up, showCategory=30) + ggtitle("GO BP Higher in MUT when REST is silenced") + 
  theme(plot.title = element_text(size = 20, face = "bold",hjust=0.9))

#for(i in 1:length(kk_nfatc1_up@result$geneID)) {
#  info <- getBM(strsplit(kk_nfatc1_up@result$geneID[i], split = "/"),  filters = "entrezgene", mart = mart, attributes = c("entrezgene", "external_gene_name"))
#  kk_nfatc1_up@result$geneID[i] <- paste(c(info$external_gene_name), collapse = "/")
#}
#o<-order(kk_nfatc1_up[,6],decreasing = FALSE)
#KEGG<-as.data.frame(kk_nfatc1_up)[o,]
#write.table(KEGG,"GO BP Higher in MUT when REST is silenced.txt", sep="\t",row.names = TRUE,col.names =NA)


GENELIST=fr3$log2FoldChange-fr4$log2FoldChange
names(GENELIST)<-res$entrezgene

w<-which(fr1$padj<0.05)
edox <- setReadable(kk_nfatc1_up, 'org.Hs.eg.db', 'ENTREZID')
###Wersja "gałązka"
cnetplot(edox, foldChange=GENELIST)
cnetplot(edox, foldChange=GENELIST) +
scale_color_gradient(low ="#adff2f" , high ="#006400" ,name="custom")

