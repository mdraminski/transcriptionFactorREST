## ---------------------------

## Purpose of script: patients survival analysis using RNA-seq data (TCGA, glioma)
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------
setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/ADria_TF_story")
meth<-readRDS("MEDIAN_methylations_TCGA_LGG_GBM_per_transcript_1kb_up_down.RDS")

setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_MS_1st_2022")
pos<-readRDS("MotifsPeaksGenes_forBartek.rds")
trans<-sapply(1:length(unique(pos$transcriptId)),function(i) strsplit(unique(pos$transcriptId),".",fixed=TRUE)[[i]][1])
m<-sapply(1:length(trans), function(i) match(trans[i],meth$ensembl_transcript_id))
##many transcripts missing
ww<-which(is.na(m))

M<-sapply(1:length(trans), function(i) match(unique(pos$ENSEMBL)[i],meth$ensembl_gene_id))

KL<-m
KL[which(is.na(m))]=M[which(is.na(m))]




meth2<-meth[na.omit(KL),14:ncol(meth)]
rownames(meth2)<-unique(pos$SYMBOL)[which(!is.na(KL))]
setwd("/media/bartek/OS/TCGA_ATRX_for_MCFS/")
clinical<-read.table("Ceccarelli_2016_clinical_table.csv", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
a<-which(clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(meth2),1,12))
atrx<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(meth2),1,12))
codel<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(meth2),1,12))
wt<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(meth2),1,12))
idh_mut<-na.omit(at)
a<-which(clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G4")
at<-match(clinical[a,1],substr(colnames(meth2),1,12))
g4_wt<-na.omit(at)





library("RColorBrewer")
library(gplots)
library(ggplot2)
##color-blinded friendly
#install.packages("viridis")
library("viridis")
POP=brewer.pal(9,"Set1")


Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)





lim<-meth2[,c(atrx,codel,idh_mut,wt,g4_wt)]
colnames(lim)<-NULL




par(cex.main=1.2) # adjust font size of titles
asd<-heatmap.2(as.matrix(lim), main = "TCGA methylation REST targets",
               keysize=0.95, #change size pf color key
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               dendrogram = "row", # no dendrogram for rows
               #Rowv = "NA", # * use self-made dendrogram
               Colv = "NA", # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               
               trace="none", # hide trace
               density.info="none", # hide histogram
               
               margins = c(0,9), # margin on top(bottom) and left(right) side.
               cexRow=1.2, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               ColSideColors = c(rep(POP[3],length(atrx)),rep(POP[5],length(codel)),rep(POP[8],length(idh_mut)),rep( "red",length(wt)),rep("#660000",length(g4_wt))) #to add nice colored strips
)

legend('left', xpd=TRUE,legend=c("ATRX-IDH","1p19q","IDH","WT","G4_WT"), col=c(POP[c(3,5,8)],"red","#660000"),pch=15,cex=0.95,bty = "n",inset=c(-0.05,0))



lim<-meth2[,c(atrx,codel,idh_mut,wt,g4_wt)]
colnames(lim)<-NULL



par(cex.main=1.2) # adjust font size of titles
asd<-heatmap.2(as.matrix(lim), main = "TCGA methylation REST targets",
               keysize=0.95, #change size pf color key
               # reorderfun=function(d, w) reorder(d, w, agglo.FUN = mean),
               # order by branch mean so the deepest color is at the top
               dendrogram = "row", # no dendrogram for rows
               #Rowv = "NA", # * use self-made dendrogram
               Colv = "NA", # make sure the columns follow data's order
               col = Colors,# color pattern of the heatmap
               
               trace="none", # hide trace
               density.info="none", # hide histogram
               
               margins = c(0,9), # margin on top(bottom) and left(right) side.
               cexRow=1.2, cexCol = 1.2, # size of row / column labels
               #xlab = "Subtype",
               srtCol=0, adjCol = c(0.5,1), # adjust the direction of row label to be horizontal
               # margin for the color key
               # ("bottom.margin", "left.margin", "top.margin", "left.margin" )
               key.par=list(mar=c(5,1,3,1)),
               ColSideColors = c(rep("purple",length(c(atrx,idh_mut,codel))),rep("red",length(wt)),rep("#660000",length(g4_wt))) #to add nice colored strips
)

legend('left', xpd=TRUE,legend=c("GII/GIII-IDHmut","GII/GIII-IDHwt","GIV"), col=c("purple","red","#660000"),pch=15,cex=0.95,bty = "n",inset=c(-0.09,0))







load("/media/bartek/OS/TCGA_GBM_LGG_DESQ_RNASEQ/Glioma_TCGA_Merged.RData")


sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
w<-sapply(1:length(nm), function(i) length(nm[[i]]))
W<-which(w<2)

data2<-data[,W]
colnames(data2)<-unique(substr(unlist(nm[W]),1,16))
###checking if there is no normal samples
WW<-which(substr(colnames(data2),14,15)<10)
data3<-data2[,WW]
data4<-log2(data3+1)

###Annotating 
library(biomaRt)
#mart<-useEnsembl(dataset="hsapiens_gene_ensembl",biomart="ensembl")
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='may2021.archive.ensembl.org')
gen<- getBM(values = rownames(data3), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description","entrezgene_id"))

m<-sapply(1:length(rownames(meth2)), function(i) which(rownames(meth2)[[i]]==gen[,2]))
M<-sapply(1:length(m), function(i) m[[i]][1])
lim_E<-data4[gen[na.omit(M),1],]

tmp<-substr(colnames(meth2),1,16)
lim_M<-meth2
colnames(lim_M)<-tmp

MM<-intersect(colnames(lim_M),colnames(lim_E))
m1<-match(MM,colnames(lim_E))
m2<-match(MM,colnames(lim_M))

LIM_E<-lim_E[,m1]
LIM_M<-lim_M[which(!is.na(M)),m2]

setwd("/media/bartek/OS/TCGA_ATRX_for_MCFS/")
clinical<-read.table("Ceccarelli_2016_clinical_table.csv", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
a<-which(clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="Mutant" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(LIM_M),1,12))
atrx<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="codel" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(LIM_M),1,12))
codel<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(LIM_M),1,12))
wt<-na.omit(at)
a<-which(clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G2" |
           clinical[,"ATRX status"]=="WT" & clinical[,"1p/19q codeletion"]=="non-codel" & clinical[,"IDH status"]=="Mutant" & clinical[,"Grade"]=="G3")
at<-match(clinical[a,1],substr(colnames(LIM_M),1,12))
idh_mut<-na.omit(at)
a<-which(clinical[,"IDH status"]=="WT" & clinical[,"Grade"]=="G4")
at<-match(clinical[a,1],substr(colnames(LIM_M),1,12))
g4_wt<-na.omit(at)


cor<-sapply(1:nrow(LIM_E),function(i) cor(as.numeric(LIM_E[i,]),as.numeric(LIM_M[i,])))
cor_idh<-sapply(1:nrow(LIM_E),function(i) cor(as.numeric(LIM_E[i,c(idh_mut,atrx,codel)]),as.numeric(LIM_M[i,c(idh_mut,atrx,codel)])))
cor_wt<-sapply(1:nrow(LIM_E),function(i) cor(as.numeric(LIM_E[i,wt]),as.numeric(LIM_M[i,wt])))
cor_g4<-sapply(1:nrow(LIM_E),function(i) cor(as.numeric(LIM_E[i,g4_wt]),as.numeric(LIM_M[i,g4_wt])))


CORR<-cbind(as.numeric(cor),as.numeric(cor_wt),as.numeric(cor_idh),as.numeric(cor_g4))
rownames(CORR)<-rownames(LIM_M)
colnames(CORR)<-c("All","GII/GIII-WT","GII/GIII-MUT","GIV")
library(corrplot)
scalebluered <- colorRampPalette(c("blue","white","orange"))(200)
corrplot(CORR,cl.pos="n",mar=c(0,0,0,3),tl.col = "black",col=scalebluered)
colorlegend(xlim=c(5,10), ylim=c(5,10), scalebluered, c(seq(-1,1,.4)), align="l", vertical=TRUE, addlabels=TRUE)


p<-sapply(1:nrow(LIM_E),function(i) cor.test(as.numeric(LIM_E[i,]),as.numeric(LIM_M[i,]))$p.value)
p_idh<-sapply(1:nrow(LIM_E),function(i) cor.test(as.numeric(LIM_E[i,c(idh_mut,atrx,codel)]),as.numeric(LIM_M[i,c(idh_mut,atrx,codel)]))$p.value)
p_wt<-sapply(1:nrow(LIM_E),function(i) cor.test(as.numeric(LIM_E[i,wt]),as.numeric(LIM_M[i,wt]))$p.value)
p_g4<-sapply(1:nrow(LIM_E),function(i) cor.test(as.numeric(LIM_E[i,g4_wt]),as.numeric(LIM_M[i,g4_wt]))$p.value)

CORR_p<-cbind(p.adjust(p),p.adjust(p_wt),p.adjust(p_idh),p.adjust(p_g4))


corrplot(CORR,cl.pos="n",mar=c(0,0,0,3),,col=scalebluered,tl.col = "black",pch.cex = 1.5,p.mat = CORR_p, sig.level = c(0.05),insig = 'label_sig')
colorlegend(xlim=c(5,10), ylim=c(5,10), scalebluered, c(seq(-1,1,.4)), align="l", vertical=TRUE, addlabels=TRUE)



setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_MS_1st_2022")
sor<-readRDS("corrplot_sort_labels.rds")

col_activ <- rep("darkorange1",14)
col_rep <- rep("forestgreen", 29)
col_labels <- c(col_activ, col_rep)

m<-match(sor[,1],rownames(CORR))
corrplot(t(CORR[na.omit(m),]),tl.col=col_labels,cl.pos="n",mar=c(0,0,0,3),,col=scalebluered,pch.cex = 1.5,p.mat = t(CORR_p[na.omit(m),]), sig.level = c(0.05),insig = 'label_sig')

