## ---------------------------

## Purpose of script: Creating beeswarm plots of expression of siREST dependent REST ChIPseq tarets in TCGA LGG/GBM dataset
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

load("RData/Glioma_TCGA_Merged.RData")


sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

###Adjust Gene expression data
logs = log2(data+1)

#biocLite("preprocessCore")
library("preprocessCore")

a<-normalize.quantiles(logs, copy=TRUE)
colnames(a)<-colnames(data)
rownames(a)<-rownames(data)

logs_norm<-a


###Annotating 
library("biomaRt")
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='jul2018.archive.ensembl.org')
gen<- getBM(values = rownames(logs_norm), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol"))

setwd("/media/bartek/OS/TCGA_ATRX_for_MCFS/")
clinical<-read.table("Ceccarelli_2016_clinical_table.csv", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
a<-which(clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(logs_norm))
atrx<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(logs_norm))
codel<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(logs_norm))
wt<-na.omit(at)

a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(logs_norm))
idh_mut<-na.omit(at)

a<-which(clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G4")
at<-match(clinical[a,1],colnames(logs_norm))
g4_wt<-na.omit(at)

a<-which(clinical[,"Grade"]=="G2")
at<-match(clinical[a,1],colnames(logs_norm))
g2<-na.omit(at)

a<-which(clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],colnames(logs_norm))
g3<-na.omit(at)

ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")

vec1<-c("SNAP25","TRIM9","CADPS","PCLO","PPFIA3","SCAMP5","STXBP5L","CDK5R2","SYN1",
        "SYP","CHRNB2","CYB5R4","ADCY8","PTPRN","VGF","UNC13A","RTN4","CDK5R1")
vec2<-c("ENO2","NMNAT2","FOXK2","PGK1","REST","EIF2AK2","PSMB2","AGPAT5","SRSF4","NSRP1")

m1<-match(vec1,gen[,2])
m2<-match(vec2,gen[,2])

M1<-match(gen[m1,1],rownames(logs_norm))
M2<-match(gen[m2,1],rownames(logs_norm))

lim<-logs_norm[M1,c(CTRL,g2,g3,g4_wt)]




library(beeswarm)


tmp<-factor(c(rep("NB",length(CTRL)),rep("GII",length(g2)),rep("GIII",length(g3)),rep("GIV",length(g4_wt))))
cell_type<-relevel(tmp,"NB")

# Bee swarm plot by group
par(mar=c(3,5,2,0),cex.axis=1.5,mfrow=c(3,6),cex.main=2)
TT<-1
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20), 
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="log2 (FPKM+1)",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-2
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-3
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-4
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-5
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-6
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,5,2,0))
TT<-7
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="log2 (FPKM+1)",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-8
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-9
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-10
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-11
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])

par(mar=c(3,3,2,1))
TT<-12
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,5,2,0))
TT<-13
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="log2 (FPKM+1)",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-14
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-15
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-16
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-17
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])
par(mar=c(3,3,2,1))
TT<-18
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec1[TT])




lim<-logs_norm[M2,c(CTRL,g2,g3,g4_wt)]

par(mar=c(3,5,2,0),cex.axis=1.5,mfrow=c(2,5),cex.main=2)
TT<-1
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="log2 (FPKM+1)",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-2
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-3
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-4
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-5
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,5,2,0))
TT<-6
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="log2 (FPKM+1)",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-7
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.75,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-8
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-9
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])
par(mar=c(3,3,2,1))
TT<-10
beeswarm(lim[TT,]~cell_type,
         pch = c(19,20,20,20),
         col = c("#006666","#FF9933","#CC0000","#660000"), bty='n',
         xlab="",ylab="",cex.lab=2,cex=0.5,main=vec2[TT])

