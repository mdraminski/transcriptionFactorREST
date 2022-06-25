## ---------------------------

## Purpose of script: Testing expression of ECM-related genes targeted by siREST in TCGA LGG/GBM datasets
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

load("/media/bartek/OS/TCGA_GBM_LGG_DESQ_RNASEQ/Glioma_TCGA_Merged.RData")


sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

###Adjust Gene expression data
##changing RPKM<1 into 1(log0 is not calculatable)

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

setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_siRNA/2ndary_NOV _2020")
kk<-read.table("GO BP Higher in MUT when REST is silenced.txt", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "\"", comment.char = "", check.names=FALSE);

######Heatmap G4 WT
r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],g4_wt])
lim<-logs_norm[,g4_wt]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:14],lim3[,130:143])

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GIV WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(15),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1))
)




####grades II and III

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],c(codel,atrx,idh_mut)])
lim<-logs_norm[,c(codel,atrx,idh_mut)]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:36],lim3[,330:365])

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-MUT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(15),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1))
)





####grades II and III IDH-WT

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],wt])
lim<-logs_norm[,wt]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:8],lim3[,72:79])

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = inferno(15),# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1))
)






####################################
####Improved 25th of jan 2021

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

######Heatmap G4 WT
r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],g4_wt])
lim<-logs_norm[,g4_wt]
lim2<-lim[,o]

vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:14],lim3[,130:143])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:14], lim4[i,15:28], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
asd<-heatmap.2(zscore2, main = 'ECM REST GIV WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)




####grades II and III

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],c(codel,atrx,idh_mut)])
lim<-logs_norm[,c(codel,atrx,idh_mut)]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:36],lim3[,330:365])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:36], lim4[i,37:72], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd


zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-MUT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=2 ##Gwiazdki istotnosci
)





####grades II and III IDH-WT

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],wt])
lim<-logs_norm[,wt]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:8],lim3[,72:79])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:8], lim4[i,9:16], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)












####################################
####Improved 25th of jan 2021 - try just 8 left/right

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

######Heatmap G4 WT
r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],g4_wt])
lim<-logs_norm[,g4_wt]
lim2<-lim[,o]

vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:8],lim3[,136:143])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:8], lim4[i,9:16], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
asd<-heatmap.2(zscore2, main = 'ECM REST GIV WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)




####grades II and III

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],c(codel,atrx,idh_mut)])
lim<-logs_norm[,c(codel,atrx,idh_mut)]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:8],lim3[,358:365])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:8], lim4[i,9:16], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd


zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-MUT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)





####grades II and III IDH-WT

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],wt])
lim<-logs_norm[,wt]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]
lim4<-cbind(lim3[,1:8],lim3[,72:79])
dim(lim4)

boby1<-sapply(1:length(rownames(lim4)), function(i) t.test(lim4[i,1:8], lim4[i,9:16], paired=FALSE)$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim4),ncol(lim4))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim4))
fd[w]="*"
labs[,ncol(lim4)]<- fd

zscore<-sapply(1:nrow(lim4), function(i) scale(lim4[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'ECM REST GII/III IDH-WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)













####################################
####Improved 25th of jan 2021 - try just 8 left/right

Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

######Heatmap G4 WT
r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],g4_wt],decreasing = TRUE)
lim<-logs_norm[,g4_wt]
lim2<-lim[,o]

vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]

boby1<-sapply(1:length(rownames(lim3)), function(i) cor.test(lim3[i,],lim3[nrow(lim3),])$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim3),ncol(lim3))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim3))
fd[w[-length(w)]]="*"
labs[,ncol(lim3)-5]<- fd

zscore<-sapply(1:nrow(lim3), function(i) scale(lim3[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
asd<-heatmap.2(zscore2, main = 'GIV WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)




####grades II and III

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],c(codel,atrx,idh_mut)],decreasing = TRUE)
lim<-logs_norm[,c(codel,atrx,idh_mut)]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]

boby1<-sapply(1:length(rownames(lim3)), function(i) cor.test(lim3[i,],lim3[nrow(lim3),])$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim3),ncol(lim3))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim3))
fd[w[-length(w)]]="*"
labs[,ncol(lim3)-10]<- fd

zscore<-sapply(1:nrow(lim3), function(i) scale(lim3[i,],center=TRUE, scale=TRUE))

zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'GII/III IDH-MUT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)





####grades II and III IDH-WT

r1<-grep("REST",gen[,2])

o<-order(logs_norm[gen[r1,1],wt],decreasing = TRUE)
lim<-logs_norm[,wt]
lim2<-lim[,o]
vec<-c(unlist(strsplit(kk$"geneID"[1],"/")),"REST")
gh<-sapply(1:length(vec), function(i) which(gen[,2]==vec[i]))

g<-gen[gh,1]
m<-match(g,rownames(logs_norm))
lim3<-lim2[m,]

boby1<-sapply(1:length(rownames(lim3)), function(i) cor.test(lim3[i,],lim3[nrow(lim3),])$p.value)
FDR1<-p.adjust(boby1,"fdr")
labs <- matrix("",nrow(lim3),ncol(lim3))
w<-which(FDR1<0.05)
fd<-rep("",nrow(lim3))
fd[w[-length(w)]]="*"
labs[,ncol(lim3)-5]<- fd

zscore<-sapply(1:nrow(lim3), function(i) scale(lim3[i,],center=TRUE, scale=TRUE))
zscore2<-t(zscore)
rownames(zscore2)<-gen[gh,2]
colnames(zscore2)<-NULL
POP=brewer.pal(9,"Set1")
par(cex.main=1.2) # adjust font size of titles
library(viridis)
asd<-heatmap.2(zscore2, main = 'GII/III IDH-WT',
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               #dendrogram = "row", # no dendrogram for columns
               Rowv = "NA", # * use self-made dendrogram
               Colv = "NA",  # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               dendrogram = "none",
               trace="none", # hide trace
               density.info="none",  # hide histogram
               
               margins = c(0,6),  # margin on top(bottom) and left(right) side.
               cexRow=NULL, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               cellnote=labs, notecol="black", notecex=3 ##Gwiazdki istotnosci
)






