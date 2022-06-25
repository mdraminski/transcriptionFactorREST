## ---------------------------

## Purpose of script: Deafine the function of REST-activated and REST-repressed genes
##
## Author: Bartosz Wojtas
##
## Date Created: 2022-03-24
##

## ---------------------------
          

          load("/media/bartek/OS/TCGA_GBM_LGG_DESQ_RNASEQ/Glioma_TCGA_Merged.RData")
          
          nam<-sapply(1:677, function(i) names(as.data.frame(glioma.tcga[[i]]$clinical)))
          hj<-sapply(1:length(nam), function(i) length(nam[[i]]))
          gh<-which(hj>0)
          
          glioma.tcga2<-glioma.tcga[gh]
          rm(glioma.tcga)
          gc(reset=TRUE)
          
          data<-sapply(1:length(glioma.tcga2), function(i) glioma.tcga2[[i]]$mRNAseq[[1]][[1]]$FPKM)
          nm<-sapply(1:length(glioma.tcga2), function(i) names(glioma.tcga2[[i]]$mRNAseq[[1]]))
          colnames(data)<-unique(substr(unlist(nm),1,12))
              
          ###Adjust Gene expression data
          logs = log2(data+1)
          
          library("preprocessCore")
          
          a<-normalize.quantiles(logs, copy=TRUE)
          colnames(a)<-colnames(data)
          rownames(a)<-rownames(data)
          
          logs_norm<-a
          
          
          ###Annotating 
          library("biomaRt")
          mart = useMart("ensembl",dataset="hsapiens_gene_ensembl",host='jul2018.archive.ensembl.org')
            gen<- getBM(values = rownames(logs_norm), filters = "ensembl_gene_id", mart = mart, attributes = c("ensembl_gene_id","hgnc_symbol"))
          
          rt<-which(gen[,2]=="REST")
          rest<-logs_norm[gen[rt,1],]
          
          setwd("/media/bartek/3d2f2a34-476a-482c-8763-6490c503b18c/REST_MS_1st_2021/")
          piki<-read.table("PeakAnno.csv", sep="\t",fill=TRUE,header=TRUE,stringsAsFactors=FALSE,quote = "\"", comment.char = "", check.names=FALSE);
          
          g<-grep("Promoter",unique(piki$annotation))
          pp<-which(piki$annotation==unique(piki$annotation)[g[1]] | piki$annotation==unique(piki$annotation)[g[2]] | piki$annotation==unique(piki$annotation)[g[3]])
          
          targets<-piki$SYMBOL[pp]
          ens<-piki$ENSEMBL[pp]
          
          m<-match(unique(ens),rownames(logs_norm))
          lim<-logs_norm[na.omit(m),]
          
          p_corr<-sapply(1:nrow(lim), function(i) cor.test(as.numeric(rest), as.numeric(lim[i,]))$p.value)
          e_corr<-sapply(1:nrow(lim), function(i) cor.test(as.numeric(rest), as.numeric(lim[i,]))$e)
          par(mfrow=c(1,1))
            p1<-hist(e_corr,xlim=c(-1,1),ylim=c(0,200),breaks=100)
            
            sam<-sample(c(1:60486)[-na.omit(m)],nrow(lim))
            test<-sapply(1:nrow(lim), function(i) cor.test(as.numeric(rest), as.numeric(logs_norm[sam[i],]))$e)
            p2<-hist(test,xlim=c(-1,1),ylim=c(0,200),breaks=100)
          
            plot(p2, col=rgb(1,0,0,1/4), xlim=c(-1,1), xlab="correlation value",main="Overlapped")
            plot(p1,col=rgb(0,0,1,1/4), xlim=c(-1,1),add=T)
            legend("topright",c("Correlation of REST and its targets","Random REST-gene pairs"),pch=15,col=c(rgb(0,0,1,1/4),rgb(1,0,1,1/4)))
            box()
            
              
                par(mfrow=c(1,1))
                p1<-hist(p.adjust(p_corr,"bonferroni"),xlim=c(0,0.05),ylim=c(0,1200),breaks=1000)
                
                sam<-sample(c(1:60486)[-na.omit(m)],nrow(lim))
                test<-sapply(1:nrow(lim), function(i) cor.test(as.numeric(rest), as.numeric(logs_norm[sam[i],]))$p.value)
                p2<-hist(p.adjust(test,"bonferroni"),xlim=c(0,0.05),breaks=1000,ylim=c(0,1200))
                
                plot(p2, col=rgb(1,0,0,1/4), xlim=c(0,0.05),ylim=c(0,1200), xlab="correlation p-value",main="Overlapped")
                plot(p1,col=rgb(0,0,1,1/4), xlim=c(0,0.05),ylim=c(0,1200),add=T)
                abline(v=0.001,col="red",lwd=2)
                text(x=0.027,y=400,"bonferroni_p_vale_of_0.001",cex=3,col="red")
                legend("topright",c("Corr p-value of REST and its targets","Random REST-gene pairs"),pch=15,col=c(rgb(0,0,1,1/4),rgb(1,0,1,1/4)))
                box()
                
          w<-which(p_corr<0.001)
          hist(e_corr[w],breaks = 100)      
            
          PL<-which(e_corr>0.15)
          MIN<-which(e_corr< -0.15)
            
              POS<-rownames(lim)[PL]
              NEG<-rownames(lim)[MIN]
            
              ##Bad approac
            m<-match(POS,piki$ENSEMBL)
            ##correct one
            g<-sapply(1:nrow(piki), function(i) which(piki$ENSEMBL[i]==POS))
            fg<-sapply(1:length(g), function(i) length(g[[i]]))
            G<-which(fg>0)
            #write.table(cbind(piki[m,],as.integer(sapply(1:nrow(piki[m,]), function(i) median(as.integer(piki[m[i],c(2:3)]))))),"REST_peaks_REST_activator.txt",sep="\t")
            write.table(cbind(piki[G,],as.integer(sapply(1:nrow(piki[G,]), function(i) median(as.integer(piki[G[i],c(2:3)]))))),"REST_peaks_REST_activator_CORRECT.txt",sep="\t")
              ####In the excel select just "Promoters" for the analysis!!!! 
            
            #bad approach    
            #m<-match(NEG,piki$ENSEMBL)
             # write.table(cbind(piki[m,],as.integer(sapply(1:nrow(piki[m,]), function(i) median(as.integer(piki[m[i],c(2:3)]))))),"REST_peaks_REST_repressor.txt",sep="\t")
            ##correct one
            g<-sapply(1:nrow(piki), function(i) which(piki$ENSEMBL[i]==NEG))
            fg<-sapply(1:length(g), function(i) length(g[[i]]))
            G<-which(fg>0)
            #write.table(cbind(piki[m,],as.integer(sapply(1:nrow(piki[m,]), function(i) median(as.integer(piki[m[i],c(2:3)]))))),"REST_peaks_REST_activator.txt",sep="\t")
            write.table(cbind(piki[G,],as.integer(sapply(1:nrow(piki[G,]), function(i) median(as.integer(piki[G[i],c(2:3)]))))),"REST_peaks_REST_repressor_CORRECT.txt",sep="\t")
            ####In the excel select just "Promoters" for the analysis!!!!   
  
  
              
              library("org.Hs.eg.db")
                  library("clusterProfiler")
              library("ggplot2")
              
              
              ###REPRESSED REST TARGETS
                eg = bitr(POS, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                ep = bitr(rownames(logs_norm), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                ego <- enrichGO(gene          = eg[,2],
                                universe      = ep[,2],
                                OrgDb         = org.Hs.eg.db,
                                ont           = "BP",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.05,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)
                barplot(ego, showCategory=8) + ggtitle("REST-activated targets")
  
              o<-order(ego[,6],decreasing = FALSE)
              KEGG<-as.data.frame(ego)[o,]
              write.table(KEGG,"GO_BP_activated_genes.txt",sep="\t")
              
              
                      ###REPRESSED REST TARGETS
                            eg = bitr(NEG, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                          ep = bitr(rownames(logs_norm), fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
                          ego <- enrichGO(gene          = eg[,2],
                                          universe      = ep[,2],
                                          OrgDb         = org.Hs.eg.db,
                                          ont           = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff  = 0.05,
                                          qvalueCutoff  = 0.05,
                                          readable      = TRUE)
                          barplot(ego, showCategory=8) + ggtitle("REST-repressed targets")
                          o<-order(ego[,6],decreasing = FALSE)
                          KEGG<-as.data.frame(ego)[o,]
                          write.table(KEGG,"GO_BP_repressed_genes.txt",sep="\t")
                          
   