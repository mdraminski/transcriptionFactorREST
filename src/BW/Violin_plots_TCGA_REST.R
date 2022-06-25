
## ---------------------------

## Purpose of script: Comparison of REST gerne expression level in normal brain (NB), grade II, III and IV gliomas on data from TCGA LGG/GBM datasets
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------


##############Load mRNA data
load("/media/bartek/OS/TCGA_GBM_LGG_DESQ_RNASEQ/Glioma_TCGA_Merged.RData")

sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

###Adjust Gene expression data

#create back-up copy
logs = log2(data+1)

#biocLite("preprocessCore")
library("preprocessCore")

a<-normalize.quantiles(logs, copy=TRUE)
colnames(a)<-colnames(data)
rownames(a)<-rownames(data)

logs_norm<-a


setwd("/media/bartek/OS/ATRX/")
clinical<-read.csv("1-s2.0-S009286741501692X-mmc2.csv", sep="\t",skip=1,fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "", comment.char = "", check.names=FALSE);
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



sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")
dip<-sapply(1:length(fd), function(i) as.vector(unlist(glioma.tcga[[i]]$clinical))[26])
G2<-which(dip=="G2")
G3<-which(dip=="G3")
G4<-which(dip=="G4")

non<-setdiff(c(1:677),c(CTRL,G2,G3,G4))

###Annotating 
library("biomaRt")
mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='jul2018.archive.ensembl.org')
gen<- getBM(values = rownames(logs_norm), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","external_gene_name","chromosome_name","start_position","end_position","percentage_gene_gc_content","gene_biotype","description"))
  
  
  kl<-match("REST",gen[,2])
  
  
  library(plotly)
  
  df <- as.data.frame(cbind(c(rep("NB",length(CTRL)),rep("GII",length(G2)),rep("GIII",length(G3)),rep("GIV",length(G4))),
                            c(logs_norm[gen[kl,1],CTRL],logs_norm[gen[kl,1],G2],logs_norm[gen[kl,1],G3],logs_norm[gen[kl,1],G4])))
  colnames(df)<-c("status","expr")
  df$expr<-as.numeric(as.character(df$expr))
  
  
  f1 <- list(
    family = "Arial",
    size = 25,
    color = "black"
  )

  col=c("#FFCC99","#003366","#FF9933","#CC0000","#00CCCC","#006666","#009999","#660000")
  pal <- c(I("#006666"),I("#FF9933"),I("#CC0000"),I("#660000"))
  p <- df %>%
    plot_ly(
      x = ~status,
      y = ~expr ,
      split = ~status,
      type = 'violin',
      color= ~status, colors = pal,
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      )
    ) %>% 
    layout(titlefont=f1,
      xaxis = list(title="",
        tickfont = f1,
        categoryorder = "array",
        categoryarray = ~status
      ),
        yaxis = list(title="log2 expr +1",
          tickfont = f1,
          titlefont=list(size=30),
        zeroline = F,
        tickfont = f1)
    )
  p
  
  oneway.test(df[,2]~df[,1])$p.value

  
  df <- as.data.frame(cbind(c(rep("WT",length(wt)),rep("IDH",length(c(idh_mut,atrx,codel)))),
                              c(logs_norm[gen[kl,1],wt],logs_norm[gen[kl,1],c(idh_mut,atrx,codel)])))
  colnames(df)<-c("status","expr")
  df$expr<-as.numeric(as.character(df$expr))
  pal <- c(I("red"),I("purple"))
  p <- df %>%
    plot_ly(
      x = ~status,
      y = ~expr ,
      split = ~status,
      type = 'violin',
      color= ~status, colors = pal,
      box = list(
        visible = T
      ),
      meanline = list(
        visible = T
      )
    ) %>% 
    layout(titlefont=f1,
           xaxis = list(title="",
                        tickfont = f1,
                        categoryorder = "array",
                        categoryarray = ~status
           ),
           yaxis = list(title="log2 expr +1",
                        tickfont = f1,
                        titlefont=list(size=30),
                        zeroline = F,
                        tickfont = f1)
    )
  p
  
  oneway.test(df[,2]~df[,1])$p.value
  t.test(df[,2]~df[,1])$p.value
  
  
  
  
  
  