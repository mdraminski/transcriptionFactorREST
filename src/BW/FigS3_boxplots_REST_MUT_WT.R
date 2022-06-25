
## ---------------------------

## Purpose of script: Creating boxplots of REST gene expression in U87 IDH-MUT vs IDH-WT; in U87 IDH-MUT vs IDH-WT with siCTRL and in U87 IDH-MUT vs IDH-WT in siREST
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


df1 <- as.data.frame(counts(dds1, normalized=TRUE))



setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_ATRX_RNAseq/")
counts<-read.table("feature_counts_REST_ATRX", sep="\t",fill=TRUE,header=TRUE,skip=1,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
mycounts<-counts[,c(7:ncol(counts))]
rownames(mycounts)<-gsub("gene:","",counts[,1])
a<-gsub("_Aligned.sortedByCoord.out.bam","",colnames(mycounts))
colnames(mycounts)<-a


###Annotating 
library("biomaRt")
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='jul2018.archive.ensembl.org')
gen<- getBM(values = rownames(mycounts), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description","entrezgene"))

WT<-grep("F11",colnames(mycounts))
MUT<-grep("F31",colnames(mycounts))

tmp<-cbind(mycounts[,WT],mycounts[,MUT])
mycounts<-tmp

library("DESeq2")
condition <- factor(c(rep("WT",3),rep("MUT",3)))
Patient.ID <- factor(c(1,2,3,1,2,3))


#coldata <- data.frame(row.names=colnames(mycounts), condition)
coldata <- data.frame(row.names=colnames(mycounts), condition, Patient.ID)

#create a DESeq object (New:multifactoral design including the patients ID to account for differences between the samples
#Fit an individual base line for each patient, so that patient-to-patient differences in expression before treatment are absorbed and you only look at the differences the condition)              
dds1 <- DESeqDataSetFromMatrix(countData=mycounts, colData=coldata, design=~Patient.ID + condition)
dds1
dds1 <- DESeq(dds1)

###get matrix of normalized counts
df2 <- as.data.frame(counts(dds1, normalized=TRUE))


m1<-match("REST",gen[,2])

log2(df[gen[m1,1],]+1)
WT_siCTRL<-intersect(grep("F11",colnames(df)),grep("_V11",colnames(df)))
WT_siREST<-intersect(grep("F11",colnames(df)),grep("_V21",colnames(df)))
MUT_siCTRL<-intersect(grep("F31",colnames(df)),grep("_V11",colnames(df)))
MUT_siREST<-intersect(grep("F31",colnames(df)),grep("_V21",colnames(df)))


par(mar=c(9,6,1,0))
m=matrix(NA,nrow=3,ncol=6)
m[1:3,1]=as.numeric(log2(df2[gen[m1,1],]+1)[WT])
m[1:3,2]=as.numeric(log2(df2[gen[m1,1],]+1)[MUT])
m[1:3,3]=as.numeric(log2(df[gen[m1,1],]+1)[WT_siCTRL])
m[1:3,4]=as.numeric(log2(df[gen[m1,1],]+1)[MUT_siCTRL])
m[1:3,5]=as.numeric(log2(df[gen[m1,1],]+1)[WT_siREST])
m[1:3,6]=as.numeric(log2(df[gen[m1,1],]+1)[MUT_siREST])
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]),as.numeric(m[,5]),as.numeric(m[,6]))
par(cex.axis=1.5,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
axis(side = 1, at = c(0,1,2,3,4,5,6),labels=c("","WT","MUT","WT_siCTRL","MUT_siCTRL","WT_siREST","MUT_siREST"),mgp=c(1, 1, 0),las=2,cex=1.5)
l<-mtext("REST log2 (FPKM+1)", side=2, line=3,cex=1.5)
stripchart(list(list[[1]],list[[2]]),vertical=TRUE,method="jitter",pch=19,col=c("red"),add=TRUE,cex=1.4)
stripchart(list(list[[NA]],list[[NA]],list[[3]],list[[4]]),vertical=TRUE,method="jitter",pch=19,col=c("gold"),add=TRUE,cex=1.4)
stripchart(list(list[[NA]],list[[NA]],list[[NA]],list[[NA]],list[[5]],list[[6]]),vertical=TRUE,method="jitter",pch=19,col=c("purple"),add=TRUE,cex=1.4)

