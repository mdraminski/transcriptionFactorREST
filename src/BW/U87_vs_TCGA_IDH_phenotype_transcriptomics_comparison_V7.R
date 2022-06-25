

## ---------------------------

## Purpose of script: Comparison of IDH-MUT and IDH-WT in TCGA gliomas and U87 cell model - comparison of up- and down-regulated genes
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------


##############Load mRNA data
load("/media/bartek/OS/TCGA_GBM_LGG_DESQ_RNASEQ/Glioma_TCGA_Merged.RData")



sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$HTseq)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
NM<-sapply(1:length(nm), function(i) nm[[i]][1])
colnames(data)<-substr(NM,1,12)



setwd("/media/bartek/OS/TCGA_ATRX_for_MCFS/")
clinical<-read.table("Ceccarelli_2016_clinical_table.csv", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
a<-which(clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(data))
mut<-na.omit(at)

a<-which(clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G2" |
           clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(data))
wt<-na.omit(at)


sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")
dip<-sapply(1:length(fd), function(i) as.vector(unlist(glioma.tcga[[i]]$clinical))[26])
G4<-which(dip=="G4")

non<-setdiff(c(1:677),c(CTRL,mut,wt,G4))

library("DESeq2")
tmp<-rep(NA,677)
tmp[CTRL]="ctrl"
tmp[mut]="mut"
tmp[wt]="wt"
tmp[G4]="G4"
tmp[non]<-"non"
condition <- factor(tmp)


#coldata <- data.frame(row.names=colnames(mycounts), condition)
coldata <- data.frame(row.names=colnames(data), condition)

#create a DESeq object (New:multifactoral design including the patients ID to account for differences between the samples
#Fit an individual base line for each patient, so that patient-to-patient differences in expression before treatment are absorbed and you only look at the differences the condition)              
dds1 <- DESeqDataSetFromMatrix(countData=data, colData=coldata, design=~condition)
dds1
dds1 <- DESeq(dds1)
fr1<-results(dds1,contrast = c("condition", "mut","wt"))



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
dds2 <- DESeqDataSetFromMatrix(countData=mycounts, colData=coldata, design=~Patient.ID + condition)
dds2
dds2 <- DESeq(dds2)
fr2<-results(dds2,contrast = c("condition", "MUT","WT"))

#####
#####takes a lot of time!!!
###can be loaded instead:
#load("/media/bartek/OS/TCGA_ATRX_for_MCFS/TCGA_LGG_GBM_DESeq2_own_normalization_IDH_MUT_vs_WT_DESeq2_analysis_DEC_2021.RData")


m<-match(rownames(fr2),rownames(fr1))
fr2_lim<-fr2[which(!is.na(m)),]
fr1_lim<-fr1[na.omit(m),]

w<-which(fr2_lim$padj<0.05)

fr1_lim2<-fr1_lim[w,]
fr2_lim2<-fr2_lim[w,]


o<-order(fr2_lim2$log2FoldChange)
par(mar=c(5,5,1,0),cex.axis=1.5)
plot(1:length(o),fr2_lim2$log2FoldChange[o],pch=19,col="red",cex=1,xlab="genes",ylab="log2FC",cex.lab=2,bty='n')
points(1:length(o),fr1_lim2$log2FoldChange[o],pch=19,col="grey",cex=0.5)
legend("topleft",c("U87","TCGA"),col=c("red","grey"),pch=19,cex=1.5,bty='n')
abline(lwd=2,h=0)
abline(lwd=2,v=2141)
text(1000,10,"61% concordance",cex=1.4)
text(3000,10,"38% concordance",cex=1.4)

length(which(fr1_lim2$log2FoldChange[o[1:2141]]<0))/length(which(fr2_lim2$log2FoldChange<0))
length(which(fr1_lim2$log2FoldChange[o[2142:length(o)]]>0))/length(which(fr2_lim2$log2FoldChange>0))


###sampling

sam<-sapply(1:10000, function(i) sample(fr1$log2FoldChange,2141))
sample<-sapply(1:10000, function(i) length(which(sam[,i]<0))/length(which(fr2_lim2$log2FoldChange<0)))

sam2<-sapply(1:10000, function(i) sample(fr1$log2FoldChange,1960))
sample2<-sapply(1:10000, function(i) length(which(sam2[,i]<0))/length(which(fr2_lim2$log2FoldChange>0)))


