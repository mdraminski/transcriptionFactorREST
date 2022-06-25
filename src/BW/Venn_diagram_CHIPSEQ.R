

## ---------------------------

## Purpose of script: Annotation and comparison of IDH-MUT and IDH-WT ChIPseq peaks
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------

library(GenomicFeatures)
library(GenomicRanges)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene


setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_peaks_U87_MUT_WT_APRIL_2020")
peakAnno_IDH_1 <- annotatePeak("U87_MUT_peaks_intersected_2rep.bed", tssRegion=c(-3000, 3000), 
                               TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno_WT_1 <- annotatePeak("U87_WT_peaks_intersected_2rep.bed", tssRegion=c(-3000, 3000), 
                              TxDb=txdb, annoDb="org.Hs.eg.db")

idh<-unique(as.data.frame(peakAnno_IDH_1)[,"SYMBOL"])
wt<-unique(as.data.frame(peakAnno_WT_1)[,"SYMBOL"])



venn.diagram(x = list(idh,wt), 
             category.names = c("IDH_peaks","WT_peaks"),
             filename = 'venn_diagramm_ChIPseq_REST.png',output = TRUE ,
             imagetype="png" ,
             height = 900 ,
             width = 900 ,
             resolution = 500,
             compression = "lzw",
             lwd = 4,
             lty = 'blank',
             fill = c('red','purple'),
             cex = 1,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.pos = 0,
             main.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             main = "ChIPseq_Venn")



venn.diagram(x = list(idh,wt), 
             category.names = c("",""),
             filename = 'venn_diagramm_ChIPseq_REST_no_names.png',output = TRUE ,
             imagetype="png" ,
             height = 900 ,
             width = 900 ,
             resolution = 500,
             compression = "lzw",
             lwd = 4,
             lty = 'blank',
             fill = c('red','purple'),
             cex = 1,
             fontface = "bold",
             fontfamily = "sans",
             cat.cex = 0.3,
             cat.pos = 0,
             main.cex = 0.5,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.fontfamily = "sans",
             main = "ChIPseq_Venn")

