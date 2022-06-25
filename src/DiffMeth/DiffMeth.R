## ---------------------------

## Purpose of script: discover differentially methylated DNA regions
##
## Author: Michal Draminski
##
## Date Created: 2020-03-22
##

## ---------------------------

library("data.table")
library("dplyr")
library("ggplot2")
library("stringi")
library("stringr")
library("yaml")
library("patchwork")
#BiocManager::install("mygene")
library("mygene")
#library("readxl")
#library("dunn.test")

#rm(list = ls()); cat("\014");
source("src/DiffMeth/DiffMeth.utils.R")
#prepare_inputDataDef("~/CytoMeth/results/methyl_results/", save_def_file=F)
cfg_file_name <- "config/config.dm.ActivRepGenes_U87_RESTpeaks.yml"

args <- commandArgs(TRUE)
if(length(args)>0){
  cfg_file_name <- args[1]
}

### CONFIGURATION ###
################################
config.DiffMeth <- readDiffMethConfig(file = cfg_file_name)
config.DiffMeth$precision <- 8
config.DiffMeth$top_beta_vs_exp_plots <- 30
#config.DiffMeth$brewer_palette <- "Dark2"
config.DiffMeth$brewer_palette <- "Set1"
config.DiffMeth$text_size <- 16
if(is.null(config.DiffMeth$pval_corr)){
  config.DiffMeth$pval_corr <- "fdr" #c("fdr", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none")
  config.DiffMeth$pval_max <- 0.05
}
if(is.null(config.DiffMeth$contexts)){
  config.DiffMeth$contexts <- c("allC","CG","nonCG")
}  
print(config.DiffMeth)

#config.DiffMeth$hgnc <- TRUE
#config.DiffMeth$overwrite_results <- TRUE


###################################################
#### START ####
start_time <- Sys.time()
# If you want to have a mask of gene promoters give the path to bed file of genes with the following columns
# format of the file chr, start, end, region, gene, strand
cat(paste0("Reading input regions from: ", config.DiffMeth$input_regions_path,"\n"))
inputRegions <- read_inputRegions(config.DiffMeth$input_regions_path)
print(inputRegions)

#### Add TSS promoters
if(config.DiffMeth$add_TSS_promoter){
  cat(paste0("Adding TSS promoters...\n"))
  inputRegions <- add_TSS_promoter(inputRegions, upstream=config.DiffMeth$TSS_promoter_upstream, downstream=config.DiffMeth$TSS_promoter_downstream, geneEnd=FALSE)
}
#### Add HGNC
inputRegions <- add_hgnc(inputRegions)
if(config.DiffMeth$hgnc)
  inputRegions$region <- inputRegions$gene <- inputRegions$hgnc
print(inputRegions)

#### INTERSECTION
# make the intersection of inputRegions with all input samples/files data - save the result files into config$intersect_path
# bedtools are needed
intersect_input_data(inputRegions, config.DiffMeth)
# and read all intersected files and filter them by config.DiffMeth$min_coverage
#methData <- readRDS(file.path(config.DiffMeth$results_path, "methData.rds"))
methData <- read_methData(input_folder=config.DiffMeth$intersect_data_path, config.DiffMeth)
methData_file_name <- file.path(config.DiffMeth$results_path, "methData.rds")
if(config.DiffMeth$overwrite_results | !file.exists(methData_file_name)){
  saveRDS(methData, methData_file_name)
}
methData <- readRDS(methData_file_name)
print(head(methData))
print(table(methData$group))
print(data.table(table(methData$sample_id)))

###################################################
#### REGIONS ####
#### CHI2 TEST - Compute differences between pairs of groups
chi2_regions_file <- file.path(config.DiffMeth$results_path, "chi2_regions.csv")
if(config.DiffMeth$overwrite_results | !file.exists(chi2_regions_file)){
  regionsInfo <- get_regions(methData)
  chi2_regions <- calcChi2MethDiff(methData, config.DiffMeth$contexts, config.DiffMeth)
  fwrite(chi2_regions, file = chi2_regions_file)
}else{
  regionsInfo <- get_regions(methData)
  cat(paste0("Loading chi2_regions file: ", chi2_regions_file),"\n")
  chi2_regions <- fread(chi2_regions_file)
}
print(chi2_regions[order(p.corrected)])
print(paste0("chi2_regions #rows: ", nrow(chi2_regions)))
print(paste0("chi2_regions significant #rows: ", nrow(chi2_regions[chi2_regions$significant==T,])))

#### BOXPLOTS & KW - Compute differences between pairs of groups
if(nrow(chi2_regions!=0)){
  save_plot <- config.DiffMeth$overwrite_results
  #### HEATMAP of p.corrected values
  gg <- plot_heatmap_pCorrected(chi2_regions, contexts=config.DiffMeth$contexts, region_type="Regions", config.DiffMeth, save_plot)
  # select only significant regions
  chi2_sig_regions <- chi2_regions[chi2_regions$significant==T,]
  #### KW & DESCRIPTIVE STATS
  groupsRegionsStats <- plot_boxplots_kw(methData, regionsInfo, chi2_sig_regions, region_type="Regions", config.DiffMeth, save_plot)
  if(config.DiffMeth$overwrite_results){
    fwrite(groupsRegionsStats$kw_results, file.path(config.DiffMeth$results_path, "groupsRegionsKWStats.csv"))
    fwrite(groupsRegionsStats$desc_stats, file.path(config.DiffMeth$results_path, "groupsRegionsDescStats.csv"))
  }
}

###################################################
#### RNA EXPRESSION #### [Optional]
if(config.DiffMeth$use_gene_expression_data){
  #### REGIONS
  chi2_sig_regions <- chi2_regions[chi2_regions$significant==T,]
  if(endsWith(tolower(config.DiffMeth$input_expression_data_path),".rds")){
    gene_expression_data <- readRDS(config.DiffMeth$input_expression_data_path)
  }else{
    gene_expression_data <- fread(config.DiffMeth$input_expression_data_path) %>% data.frame()
  }
  #if input data inputRegions$gene are HGNC and expression is not (ENSG) then convert expression to HGNC
  if(!startsWith(toupper(inputRegions$gene[1]),"ENSG") & startsWith(toupper(names(gene_expression_data)[2]),"ENSG")){
    expr_cols <- intersect(inputRegions$ensembl, names(gene_expression_data))
    gene_expression_data <- setDT(gene_expression_data)[,c("sample_id",expr_cols),with=FALSE]
    inputRegions <- as.data.frame(inputRegions)
    rownames(inputRegions) <- inputRegions$ensembl
    colnames(gene_expression_data) <- c("sample_id",inputRegions[expr_cols,]$hgnc)
    rownames(inputRegions) <- NULL
  }  
  print(gene_expression_data[1:6,1:6])
  print(dim(gene_expression_data))
  
  #meanBeta_regions <- calc_mean_beta(methData, chi2_regions[!is.na(chi2_regions$pval),], gene_expression_data)
  meanBeta_regions <- calc_mean_beta(methData, chi2_regions, gene_expression_data)
  if(config.DiffMeth$overwrite_results){
    saveRDS(meanBeta_regions, file.path(config.DiffMeth$results_path, "meanBeta_regions.rds"))
  }
  save_plot <- config.DiffMeth$overwrite_results
  
  # filter to the only significant regions in specific context
  meanBeta_regions <- meanBeta_regions[paste0(meanBeta_regions$context,"_",meanBeta_regions$region) %in% unique(paste0(chi2_sig_regions$context,"_",chi2_sig_regions$region))]
  #### PLOT HEATMAPS for all contexts and significant regions
  ret <- lapply(split(meanBeta_regions, by="context", drop=TRUE), 
                function(x){plot_region_sample_heatmap(x, "Regions", config.DiffMeth, save_plot, plot_size=12, plot_margin=6, plot_gene=F)})
  #### CALC slopeStatsRegions
  statsRegions <- calcSlopeStats(meanBeta_regions, methData)
  #if select all sig_regions
  #selected_region_context <- apply(unique(chi2_sig_regions[,c("region","context")]), 1, paste, collapse=".")
  selected_region_context <- unique(apply(head(statsRegions$slope[,c("region","context")],config.DiffMeth$top_beta_vs_exp_plots), 1, paste, collapse="."),
         apply(statsRegions$corr[statsRegions$corr$corr_pval<config.DiffMeth$max_pval,c("region","context")], 1, paste, collapse="."))
  #### PLOT BETA VS EXPR - general and detailed
  meanBeta_regions <- split(meanBeta_regions,list(meanBeta_regions$region, meanBeta_regions$context), drop=TRUE)
  ret <- lapply(meanBeta_regions[selected_region_context], function(x){plot_beta_vs_exp(x, "Regions", config.DiffMeth, save_plot);return(TRUE)})
  statsRegions <- statsRegions_rollup(statsRegions, selected_region_context)
  meanBeta_regions <- rbindlist(meanBeta_regions)
  if(config.DiffMeth$overwrite_results){
    fwrite(statsRegions, file.path(config.DiffMeth$results_path, "slopeStatsRegions.csv"))
  }
  
}
###################################################
###################################################

stop_time <- Sys.time()
cat(paste0("DiffMeth calculations are finished. Total time: ", format(stop_time - start_time, digits=3) ,"\n"))

