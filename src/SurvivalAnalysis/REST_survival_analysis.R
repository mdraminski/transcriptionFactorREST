## ---------------------------

## Purpose of script: patients survival analysis using RNA-seq data (TCGA, glioma)
##
## Author: Adria-Jaume Roura
##
## Date Created: 2019-10-31
##

## ---------------------------

library("preprocessCore")
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(dplyr)
## Loading TCGA (glioma dataset) clinical information and gene expression 
load("RData/Glioma_TCGA_Merged.RData")

######### 1. TCGA Data preparation #########
## Select glioma samples where we have RNAseq data
sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)
data<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$FPKM) #obtaining the FPKMs (paired-end seq)
nm<-sapply(1:length(fd), function(i) names(glioma.tcga[[i]]$mRNAseq[[1]]))
colnames(data)<-unique(substr(unlist(nm),1,12))

## Adjust Gene expression data, changing RPKM<1 into 1 (log0 is not calculatable)
a = (data < 1)
tuned_RPKMs = data
tuned_RPKMs[a] = 1
logs = log2(tuned_RPKMs)
a<-normalize.quantiles(logs, copy=TRUE)
colnames(a)<-colnames(data)
rownames(a)<-rownames(data)
logs_norm <- a
sd<-sapply(1:length(glioma.tcga), function(i) length(glioma.tcga[[i]]$mRNAseq))
fd<-which(sd>0)

## Normal brian samples and glioma grades
ctrl<-sapply(1:length(fd), function(i) glioma.tcga[[i]]$mRNAseq[[1]][[1]]$isControl)
CTRL<-which(ctrl=="YES")
dip<-sapply(1:length(fd), function(i) as.vector(unlist(glioma.tcga[[i]]$clinical))[26]) #column containing the tumour's grade 
G2<-which(dip=="G2")
G3<-which(dip=="G3")
G4<-which(dip=="G4")
HGG<-which(dip == "G3" | dip == "G4")
all_grades <- which(dip == "G2" | dip == "G3" | dip == "G4")
G2G3 <- which(dip == "G2" | dip == "G3")
non<-setdiff(c(1:677),c(CTRL,G2,G3,G4)) #not defined grades

######### 2. mRNA expression levels of REST #############
mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host = "jan2019.archive.ensembl.org")
genes <- getBM(values = c("REST"),  filters = "hgnc_symbol", mart = mart, attributes = c("ensembl_gene_id")) 
logs_norm_REST <- logs_norm[genes$ensembl_gene_id,];logs_norm_REST <- as.data.frame(logs_norm_REST)
logs_norm_REST$group <- "group"

## Rename
logs_norm_REST$group[non] <- "non"
logs_norm_REST$group[CTRL] <- "NB"
logs_norm_REST$group[G2] <- "GII"
logs_norm_REST$group[G3] <- "GIII"
logs_norm_REST$group[G4] <- "GIV"

######### 3. Survival analysis #############
mRNA_REST <- data[genes$ensembl_gene_id,]
mRNA_REST <- mRNA_REST[G2G3] #option1 focusing on GII/GIII patients: Figure 1D
mRNA_REST <- mRNA_REST[G4] #option2 focusing on GIV patients: Figure 1C
med <- median(mRNA_REST) #compute the median expression of REST
mRNA_REST <- as.data.frame(mRNA_REST);mRNA_REST$Case <- rownames(mRNA_REST)

## Read TCGA clinical data (Verhaak et al.) to include the IDH status per patient
REST_clinical <- read.table(paste(wd,"TableS1.PatientData.20151020.v3.csv", sep = ""),header = TRUE, sep = "\t")
REST_clinical <- REST_clinical[,c(1,14,17,18,22)] #select patient ID, Grade, Survival (months), Vital status (0, 1) and IDH status (WT, mut)
final <- merge(mRNA_REST, REST_clinical, by="Case") #merge RNA-seq and clinical data
final <- final[complete.cases(final), ] #Remove patients with NA information

## DF to save
#colnames(final) <- c("Case", "mRNA_REST", "Grade", "Survival_months", "Vital_status", "IDH.status")
#final$group <- ""
#final$group <- ifelse(test = final$mRNA_REST < 1.553623, "low_REST_group", no = "not_selected")
#final$group <- ifelse(test = final$mRNA_REST > 3.236715, "high_REST_group", no = final$group)
#write.table(final, file = "/media/adria/803c004b-f1cf-461f-8e10-ca16ad422ae3/REST_project/TCGA_REST_clinical_information.bed", quote = FALSE, sep = "\t",row.names = FALSE, col.names = TRUE)

## Plot REST differences based on IDH status (data not published)
p <- ggboxplot(final, x = "IDH.status", y = "mRNA_REST", 
               color = "IDH.status",
               palette = "Set1", add = "jitter",
               order = c("WT", "Mutant"),
               ylab = "log2(FPKM)", xlab = "",
               title = "TCGA (G2/G3) REST expression in IDH1wt and IDH1mut",
               legend = "right",
               #caption = "P-value based on Wilcoxon signed-rank test",
               outline = TRUE,
               outlier.shape =NA)
my_comparisons <- list(c("WT", "Mutant"))

p + stat_compare_means( comparisons = my_comparisons, method = "wilcox.test") +
  font("title", size = 20, color = "black", face = "bold") +
  font("xlab", size = 14, color = "black", face = "bold") +
  font("ylab", size = 14, color = "black", face = "bold") +
  font("caption", size = 9, color = "black", face = "bold.italic") + border()

## 1) Calculate KM on GII/GIII with IDH mutant phenotype
final <- final[final[,6] == "Mutant",]
rownames(final) <- final$Case
final[,c(3,4,5,6)] <- NULL
final$Case <- NULL
final <- t(final)
mRNA_REST <- final

## Stratify patient cohort based on high and low REST mRNA expression using computed median
f2<-which(mRNA_REST>1.25*med) 
f0.5<-which(mRNA_REST<0.65*med) 
mRNA_REST <- t(mRNA_REST)
mRNA_REST <- as.data.frame(mRNA_REST)
REST_bottom <- subset(mRNA_REST, mRNA_REST <= 0.65*med)
REST_top <- subset(mRNA_REST, mRNA_REST >= 1.25*med)

## Obtain futime and fustat data from TCGA
#BOTTOM PATIENTS
days=list()
grades=list()
follow_up=list()
counter = 0
for (id in rownames(REST_bottom)){
  counter = counter + 1
  days[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(6)] 
  follow_up[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(7)]
  grades[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(26)]
}
REST_bottom$days_to_death <- as.numeric(as.character(days))
REST_bottom$grade <- as.character(grades)
REST_bottom$follow_up <- as.numeric(as.character(follow_up))
REST_bottom$fustat <- "none"
REST_bottom <- REST_bottom %>% mutate(fustat = replace(fustat, is.na(days_to_death), 0)) #add censored data
REST_bottom <- REST_bottom %>% mutate(fustat = replace(fustat, !is.na(days_to_death), 1)) #add death status
REST_bottom <- REST_bottom[!REST_bottom$grade == "NULL",] #remove NULLS

## Replace NA values
select <- is.na(REST_bottom$follow_up)
REST_bottom[select,"follow_up"] <- REST_bottom[select,"days_to_death"]
REST_bottom <- REST_bottom %>% mutate(days_to_death = replace(days_to_death, is.na(days_to_death), 0))

## Modify cases where [days to death > follow_up] and treat them as uncensored patients
select <- REST_bottom$days_to_death > REST_bottom$follow_up
REST_bottom[select,"follow_up"] <- REST_bottom[select,"days_to_death"]

#TOP PATIENTS
days=list()
grades=list()
follow_up=list()
counter = 0
for (id in rownames(REST_top)){
  counter = counter + 1
  days[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(6)] 
  follow_up[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(7)]
  grades[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(26)]
}
REST_top$days_to_death <- as.numeric(as.character(days))
REST_top$grade <- as.character(grades)
REST_top$follow_up <- as.numeric(as.character(follow_up))
REST_top$fustat <- "none"
REST_top <- REST_top %>% mutate(fustat = replace(fustat, is.na(days_to_death), 0)) #add censored data
REST_top <- REST_top %>% mutate(fustat = replace(fustat, !is.na(days_to_death), 1)) #add death status
REST_top <- REST_top[!REST_top$grade == "NULL",] #remove NULLS

## Replace NA values
select <- is.na(REST_top$follow_up)
REST_top[select,"follow_up"] <- REST_top[select,"days_to_death"]
REST_top <- REST_top %>% mutate(days_to_death = replace(days_to_death, is.na(days_to_death), 0))

## Modify cases where [days to death > follow_up] and treat them as uncensored patients
select <- REST_top$days_to_death > REST_top$follow_up
REST_top[select,"follow_up"] <- REST_top[select,"days_to_death"]

## Merge
REST_top$rx <- "high"; REST_bottom$rx <- "low"
colnames(REST_top) <- c("mRNA_RPKM", "days_to_death", "grade", "futime", "fustat", "rx"); colnames(REST_bottom) <- c("mRNA_RPKM", "days_to_death", "grade", "futime", "fustat", "rx")
REST_top$fustat <- as.numeric(REST_top$fustat); REST_bottom$fustat <- as.numeric(REST_bottom$fustat)

## Survfit
gene_list <- c("REST")
glioma <- rbind(REST_bottom, REST_top); glioma$fustat <- as.numeric(glioma$fustat)
surv_object <- Surv(time = glioma$futime, event = glioma$fustat)
fit1 <- survfit(surv_object ~ rx, data = glioma)

## Plotting Kaplan-Meier plot for GII/GIII gliomas with IDH mutant phenotype (Figure 1D)
pdf(file = paste(wd,"Figure_1D.pdf", sep = ""), width = 10, height = 8)
ggsurvplot(fit1, data = glioma, pval = TRUE,palette = c("firebrick","dodgerblue3"), pval.coord = c(2500, 0.95), surv.median.line = "hv",
           xlab = "Time (days)", legend.labs=c(paste("High",gene_list[1],"mRNA"), paste("Low",gene_list[1],"mRNA")), legend.title = paste(gene_list[1], "expression"),
           risk.table = FALSE, conf.int = FALSE, tables.height = 0.2, ggtheme = theme_bw(), font.main = c(16, "bold", "darkblue"),
           font.x = c(16, "bold", "black"),font.y = c(16, "bold", "black"),font.tickslab = c(16, "plain", "black"), pval.size=5)
dev.off()


## 2) Calculate KM on GIV with IDH wild-type phenotype
final <- final[final[,6] == "WT",]
rownames(final) <- final$Case
final[,c(3,4,5,6)] <- NULL
final$Case <- NULL
final <- t(final)
mRNA_REST <- final

## Stratify patient cohort based on high and low REST mRNA expression using computed median
f2<-which(mRNA_REST>1.25*med) 
f0.5<-which(mRNA_REST<0.65*med) 
mRNA_REST <- t(mRNA_REST)
mRNA_REST <- as.data.frame(mRNA_REST)
REST_bottom <- subset(mRNA_REST, mRNA_REST <= 0.65*med)
REST_top <- subset(mRNA_REST, mRNA_REST >= 1.25*med)

## Obtain futime and fustat data from TCGA
#BOTTOM PATIENTS
days=list()
grades=list()
follow_up=list()
counter = 0
for (id in rownames(REST_bottom)){
  counter = counter + 1
  days[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(6)] 
  follow_up[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(7)]
  grades[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(26)]
}
REST_bottom$days_to_death <- as.numeric(as.character(days))
REST_bottom$grade <- as.character(grades)
REST_bottom$follow_up <- as.numeric(as.character(follow_up))
REST_bottom$fustat <- "none"
REST_bottom <- REST_bottom %>% mutate(fustat = replace(fustat, is.na(days_to_death), 0)) #add censored data
REST_bottom <- REST_bottom %>% mutate(fustat = replace(fustat, !is.na(days_to_death), 1)) #add death status
REST_bottom <- REST_bottom[!REST_bottom$grade == "NULL",] #remove NULLS

## Replace NA values
select <- is.na(REST_bottom$follow_up)
REST_bottom[select,"follow_up"] <- REST_bottom[select,"days_to_death"]
REST_bottom <- REST_bottom %>% mutate(days_to_death = replace(days_to_death, is.na(days_to_death), 0))

## Modify cases where [days to death > follow_up] and treat them as uncensored patients
select <- REST_bottom$days_to_death > REST_bottom$follow_up
REST_bottom[select,"follow_up"] <- REST_bottom[select,"days_to_death"]

#TOP PATIENTS
days=list()
grades=list()
follow_up=list()
counter = 0
for (id in rownames(REST_top)){
  counter = counter + 1
  days[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(6)] 
  follow_up[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(7)]
  grades[counter] <- as.vector(unlist(glioma.tcga[[id]]$clinical))[c(26)]
}
REST_top$days_to_death <- as.numeric(as.character(days))
REST_top$grade <- as.character(grades)
REST_top$follow_up <- as.numeric(as.character(follow_up))
REST_top$fustat <- "none"
REST_top <- REST_top %>% mutate(fustat = replace(fustat, is.na(days_to_death), 0)) #add censored data
REST_top <- REST_top %>% mutate(fustat = replace(fustat, !is.na(days_to_death), 1)) #add death status
REST_top <- REST_top[!REST_top$grade == "NULL",] #remove NULLS

## Replace NA values
select <- is.na(REST_top$follow_up)
REST_top[select,"follow_up"] <- REST_top[select,"days_to_death"]
REST_top <- REST_top %>% mutate(days_to_death = replace(days_to_death, is.na(days_to_death), 0))

## Modify cases where [days to death > follow_up] and treat them as uncensored patients
select <- REST_top$days_to_death > REST_top$follow_up
REST_top[select,"follow_up"] <- REST_top[select,"days_to_death"]

## Merge
REST_top$rx <- "high"; REST_bottom$rx <- "low"
colnames(REST_top) <- c("mRNA_RPKM", "days_to_death", "grade", "futime", "fustat", "rx"); colnames(REST_bottom) <- c("mRNA_RPKM", "days_to_death", "grade", "futime", "fustat", "rx")
REST_top$fustat <- as.numeric(REST_top$fustat); REST_bottom$fustat <- as.numeric(REST_bottom$fustat)

## Survfit
gene_list <- c("REST")
glioma <- rbind(REST_bottom, REST_top); glioma$fustat <- as.numeric(glioma$fustat)
surv_object <- Surv(time = glioma$futime, event = glioma$fustat)
fit1 <- survfit(surv_object ~ rx, data = glioma)

## Plotting Kaplan-Meier plot for GIV gliomas with IDH mutant wild-type (Figure 1C)
pdf(file = paste(wd,"Figure_1C.pdf", sep = ""), width = 10, height = 8)
ggsurvplot(fit1, data = glioma, pval = TRUE,palette = c("firebrick","dodgerblue3"), pval.coord = c(700, 0.95), surv.median.line = "hv",
           xlab = "Time (days)", legend.labs=c(paste("High",gene_list[1],"mRNA"), paste("Low",gene_list[1],"mRNA")), legend.title = paste(gene_list[1], "expression"),
           risk.table = FALSE, conf.int = FALSE, tables.height = 0.2, ggtheme = theme_bw(), font.main = c(16, "bold", "darkblue"),
           font.x = c(16, "bold", "black"),font.y = c(16, "bold", "black"),font.tickslab = c(16, "plain", "black"), pval.size=5)
dev.off()


####### END #######


