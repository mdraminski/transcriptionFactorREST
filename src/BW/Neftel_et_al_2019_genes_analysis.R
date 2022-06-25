## ---------------------------

## Purpose of script: Finding what are the correlations of selected siRNA-dependent REST ChIPseq targets in scRNAseq data from Neftel et al 2019
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/Salwador_Oncorndi/single_cell_publicDATA/")
CHI3L1<-read.table("IDHwtGBM.processed.SS2.logTPM_genes_Fig6_REST.txt", sep="\t",fill=TRUE,header=FALSE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
tmp<-CHI3L1[,1]
CHI3L1<-CHI3L1[2:length(CHI3L1)]
rownames(CHI3L1)<-tmp
names<-read.table("names_columns.txt", sep="\t",fill=TRUE,header=FALSE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
names<-names[2:length(names)]
metadata<-read.table("IDHwt.GBM.Metadata.SS2.txt", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
metadata<-metadata[-1,]

vec1<-c("SNAP25","TRIM9","CADPS","PCLO","PPFIA3","SCAMP5","STXBP5L","CDK5R2","SYN1",
        "SYP","CHRNB2","CYB5R4","ADCY8","PTPRN","VGF","UNC13A","RTN4","CDK5R1")
vec2<-c("ENO2","NMNAT2","FOXK2","PGK1","REST","EIF2AK2","PSMB2","AGPAT5","SRSF4","NSRP1")

V1<-sapply(1:length(vec1), function(i) which(rownames(CHI3L1)==vec1[i]))
V2<-sapply(1:length(vec2), function(i) which(rownames(CHI3L1)==vec2[i]))




vec<-unique(metadata$CellAssignment)
malignant<-which(metadata$CellAssignment==vec[1])
macrophage<-which(metadata$CellAssignment==vec[2])
oligo<-which(metadata$CellAssignment==vec[3])
tcell<-which(metadata$CellAssignment==vec[4])









TT<-1
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,6,4,1),mfrow=c(5,6))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("logTPM", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-2
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-3
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-4
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-5
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-6
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-7
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,6,4,1))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("logTPM", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-8
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-9
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-10
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-11
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-12
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-13
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,6,4,1))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("logTPM", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-14
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-15
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-16
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-17
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-18
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V1[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V1[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V1[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec1[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)



TT<-1
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,6,4,1))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("logTPM", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-2
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-3
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-4
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-5
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-6
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-7
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,6,4,1))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("logTPM", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)

TT<-8
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-9
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)


TT<-10
m1<-na.omit(as.numeric(CHI3L1[V2[TT],malignant]))
M1<-m1[which(m1!=0)]
m2<-na.omit(as.numeric(CHI3L1[V2[TT],macrophage]))
M2<-m2[which(m2!=0)]
m3<-na.omit(as.numeric(CHI3L1[V2[TT],oligo]))
M3<-m3[which(m3!=0)]
m4<-na.omit(as.numeric(CHI3L1[V2[TT],tcell]))
M4<-m4[which(m4!=0)]
par(mar=c(5,2,4,2))
m=matrix(NA,nrow=8000,ncol=4)
m[1:length(M1),1]=M1
m[1:length(M2),2]=M2
m[1:length(M3),3]=M3
m[1:length(M4),4]=M4
list<-list(as.numeric(m[,1]),as.numeric(m[,2]),as.numeric(m[,3]),as.numeric(m[,4]))
colnames(m)=vec
par(cex.axis=2,family="Arial")
boxplot(data.frame(m),cex.lab=2.5,cex=3,main="",outpch=NA,xaxt='n',frame="F")
title(vec2[TT],adj  = 0,cex.main=1.8)
par(srt=90)
axis(side = 1, at = c(0,1,2,3,4),labels=c("","GBM","Mp","Oligo","TCell"),mgp=c(1, 1, 0),las=2)
l<-mtext("", side=2, line=3,cex=1.5)
stripchart(list,vertical=TRUE,method="jitter",pch=19,col=c("orange"),add=TRUE,cex=0.4)





wd<-"/home/bartek/Pulpit/podsumowanie_RES_V7_ALL_FIGURES/"
pdf(paste(wd, "Neftel_correlations_V7.pdf", sep = ""), onefile = T, height = 10, width = 10)
TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"MESlike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"MESlike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"MESlike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"MESlike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"MESlike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)







TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"MESlike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"MESlike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"MESlike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"MESlike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"MESlike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)






TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"NPClike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"NPClike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"NPClike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"NPClike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"NPClike1"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike1[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)






TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"NPClike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"NPClike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"NPClike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"NPClike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"NPClike2"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$NPClike2[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)





TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"AClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"AClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"AClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"AClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"AClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$AClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)





TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"OPClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"OPClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"OPClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"OPClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"OPClike"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$OPClike[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)





TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"G1S"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"G1S"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"G1S"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"G1S"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"G1S"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G1S[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)





TT<-1
par(mar=c(5,6,1,0),mfrow=c(5,6))
a<-vec1[TT]
b<-"G2M"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)


TT<-4
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"G2M"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-11
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-12
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-13
par(mar=c(5,6,1,0))
a<-vec1[TT]
b<-"G2M"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-14
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-15
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-16
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-17
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-18
par(mar=c(5,4,1,1))
a<-vec1[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-1
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"G2M"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-2
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-3
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-4
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-5
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-6
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-7
par(mar=c(5,6,1,0))
a<-vec2[TT]
b<-"G2M"
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-8
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-9
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)

TT<-10
par(mar=c(5,4,1,1))
a<-vec2[TT]
b<-""
m1<-na.omit(as.numeric(CHI3L1[V1[TT],malignant]))
x1<-as.numeric(CHI3L1[V2[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$G2M[malignant[which(m1!=0)]])
BB<-if (cor.test(x1,y1)$p.value<0.001) {
  "darkorchid" } else { "darkgrey" }
plot(x1,y1, xlab=a, ylab=b, cex.lab=2,cex=0.5,bty='n')
ff<-lm(y1~x1)
points(y1~x1, cex=0.2,pch=20,col="orange")
abline(ff, lwd=7,col=BB)
dev.off()








V1<-sapply(1:length(vec1), function(i) which(rownames(CHI3L1)==vec1[i]))

vec<-unique(metadata$CellAssignment)
malignant<-which(metadata$CellAssignment==vec[1])
macrophage<-which(metadata$CellAssignment==vec[2])
oligo<-which(metadata$CellAssignment==vec[3])
tcell<-which(metadata$CellAssignment==vec[4])

x1<-as.numeric(CHI3L1[V1[TT],])[malignant[which(m1!=0)]]
y1<-as.numeric(metadata$MESlike1[malignant[which(m1!=0)]])

X<-sapply(1:length(V1), function(i) as.numeric(CHI3L1[V1[i],])[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]])
Y_mes1<-sapply(1:length(V1), function(i) as.numeric(metadata$MESlike1[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_mes2<-sapply(1:length(V1), function(i) as.numeric(metadata$MESlike2[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_npc1<-sapply(1:length(V1), function(i) as.numeric(metadata$NPClike1[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_npc2<-sapply(1:length(V1), function(i) as.numeric(metadata$NPClike2[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_ac<-sapply(1:length(V1), function(i) as.numeric(metadata$AClike[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_opc<-sapply(1:length(V1), function(i) as.numeric(metadata$OPClike[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_g1s<-sapply(1:length(V1), function(i) as.numeric(metadata$G1S[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_g2m<-sapply(1:length(V1), function(i) as.numeric(metadata$G2M[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))


cor_mes1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes1[[i]])$e)
p_mes1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes1[[i]])$p.value)
cor_mes2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes2[[i]])$e)
p_mes2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes2[[i]])$p.value)

cor_npc1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc1[[i]])$e)
p_npc1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc1[[i]])$p.value)
cor_npc2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc2[[i]])$e)
p_npc2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc2[[i]])$p.value)


cor_ac<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_ac[[i]])$e)
p_ac<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_ac[[i]])$p.value)
cor_opc<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_opc[[i]])$e)
p_opc<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_opc[[i]])$p.value)

cor_g1s<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g1s[[i]])$e)
p_g1s<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g1s[[i]])$p.value)
cor_g2m<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g2m[[i]])$e)
p_g2m<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g2m[[i]])$p.value)

COR<-cbind(cor_mes1,cor_mes2,cor_npc1,cor_npc2,cor_ac,cor_opc,cor_g1s,cor_g2m)
rownames(COR)<-vec1
colnames(COR)<-c("MESlike1","MESlike2","NPClike1","NPClike2","AClike","OPClike","G1S","G2M")



library(corrplot)
library(colorspace)
library(RColorBrewer)
scalebluered <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(100))

CORR_p<-cbind(p.adjust(p_mes1),p.adjust(p_mes2),p.adjust(p_npc1),p.adjust(p_npc2),p.adjust(p_ac),p.adjust(p_opc),p.adjust(p_g1s),p.adjust(p_g2m))

corrplot(t(COR),method="color",cl.pos="n",mar=c(0,0,0,3),tl.col = "black", p.mat = t(CORR_p), sig.level = c(0.05),insig = 'label_sig',
         pch.cex = 1.1,col=scalebluered)
#colorlegend(xlim=c(5,10), ylim=c(5,10), scalebluered, c(seq(-1,1,.5)), align="l", vertical=TRUE, addlabels=TRUE)

corrplot(COR,method="color",cl.pos="n",mar=c(0,0,0,3),tl.col = "black", p.mat = CORR_p, sig.level = c(0.05),insig = 'label_sig',
         pch.cex = 1.1,col=scalebluered)
#colorlegend(xlim=c(5,10), ylim=c(5,10), scalebluered, c(seq(-1,1,.5)), align="l", vertical=TRUE, addlabels=TRUE)





V1<-sapply(1:length(vec2), function(i) which(rownames(CHI3L1)==vec2[i]))

vec<-unique(metadata$CellAssignment)
malignant<-which(metadata$CellAssignment==vec[1])
macrophage<-which(metadata$CellAssignment==vec[2])
oligo<-which(metadata$CellAssignment==vec[3])
tcell<-which(metadata$CellAssignment==vec[4])


X<-sapply(1:length(V1), function(i) as.numeric(CHI3L1[V1[i],])[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]])
Y_mes1<-sapply(1:length(V1), function(i) as.numeric(metadata$MESlike1[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_mes2<-sapply(1:length(V1), function(i) as.numeric(metadata$MESlike2[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_npc1<-sapply(1:length(V1), function(i) as.numeric(metadata$NPClike1[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_npc2<-sapply(1:length(V1), function(i) as.numeric(metadata$NPClike2[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_ac<-sapply(1:length(V1), function(i) as.numeric(metadata$AClike[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_opc<-sapply(1:length(V1), function(i) as.numeric(metadata$OPClike[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_g1s<-sapply(1:length(V1), function(i) as.numeric(metadata$G1S[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))
Y_g2m<-sapply(1:length(V1), function(i) as.numeric(metadata$G2M[malignant[which(na.omit(as.numeric(CHI3L1[V1[i],malignant]))!=0)]]))


cor_mes1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes1[[i]])$e)
p_mes1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes1[[i]])$p.value)
cor_mes2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes2[[i]])$e)
p_mes2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_mes2[[i]])$p.value)

cor_npc1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc1[[i]])$e)
p_npc1<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc1[[i]])$p.value)
cor_npc2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc2[[i]])$e)
p_npc2<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_npc2[[i]])$p.value)


cor_ac<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_ac[[i]])$e)
p_ac<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_ac[[i]])$p.value)
cor_opc<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_opc[[i]])$e)
p_opc<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_opc[[i]])$p.value)

cor_g1s<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g1s[[i]])$e)
p_g1s<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g1s[[i]])$p.value)
cor_g2m<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g2m[[i]])$e)
p_g2m<-sapply(1:length(X), function(i) cor.test(X[[i]],Y_g2m[[i]])$p.value)

COR<-cbind(cor_mes1,cor_mes2,cor_npc1,cor_npc2,cor_ac,cor_opc,cor_g1s,cor_g2m)
rownames(COR)<-vec2
colnames(COR)<-c("MESlike1","MESlike2","NPClike1","NPClike2","AClike","OPClike","G1S","G2M")



library(corrplot)
library(colorspace)
library(RColorBrewer)
scalebluered <- rev(colorRampPalette(brewer.pal(8, "RdBu"))(100))

CORR_p<-cbind(p.adjust(p_mes1),p.adjust(p_mes2),p.adjust(p_npc1),p.adjust(p_npc2),p.adjust(p_ac),p.adjust(p_opc),p.adjust(p_g1s),p.adjust(p_g2m))

corrplot(t(COR),method="color",cl.pos="n",mar=c(0,0,0,3),tl.col = "black", p.mat = t(CORR_p), sig.level = c(0.05),insig = 'label_sig',
         pch.cex = 1.1,col=scalebluered)
#colorlegend(xlim=c(2.5,5), ylim=c(2.5,5), scalebluered, c(seq(-1,1,.5)), align="l", vertical=TRUE, addlabels=TRUE)

corrplot(COR,method="color",cl.pos="n",mar=c(0,0,0,3),tl.col = "black", p.mat = CORR_p, sig.level = c(0.05),insig = 'label_sig',
         pch.cex = 1.1,col=scalebluered)
