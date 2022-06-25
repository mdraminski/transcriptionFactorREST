## ---------------------------

## Purpose of script: discover differentially methylated DNA regions
##
## Author: Michal Draminski
##
## Date Created: 2020-03-22
##

## ---------------------------

#rm(list = ls()); cat("\014");
source("src/DiffMeth/DiffMethCytoRegions.utils.R")

###################################################
#### CYTO REGIONS #### [Optional]
if(config.DiffMeth$cyto_regions_find){
  chi2_cytoRegions_file <- file.path(config.DiffMeth$results_path, "chi2_results_cytoRegions.csv")
  methData_cytoRegions_file <- file.path(config.DiffMeth$results_path, "methData_cytoRegions.rds")
  #### CHI2 TEST - Compute differences between pairs of groups
  if(config.DiffMeth$overwrite_results | !file.exists(chi2_cytoRegions_file)){
    cytoRegionsInfo <- find_cytoRegions(methData, config.DiffMeth$cyto_regions_max_distance)
    methData_cytoRegions <- cytoRegions_merge_methData(cytoRegionsInfo, methData, config.DiffMeth)
    check_methData_cytoRegions(methData_cytoRegions)
    saveRDS(methData_cytoRegions, methData_cytoRegions_file)
    # Compute differences between pairs of groups of cyto_regions
    chi2_cytoRegions <- calcChi2MethDiff(methData_cytoRegions, config.DiffMeth$contexts, config.DiffMeth)
    fwrite(chi2_cytoRegions, file = chi2_cytoRegions_file)
  }else{
    cat(paste0("Loading methData_cytoRegions from file: ", methData_cytoRegions_file),"\n")
    methData_cytoRegions <- readRDS(methData_cytoRegions_file)
    print(methData_cytoRegions)
    cytoRegionsInfo <- get_cytoRegions(methData_cytoRegions)
    cat(paste0("Loading chi2_cyto_regions file: ", chi2_cytoRegions_file),"\n")
    chi2_cytoRegions <- fread(chi2_cytoRegions_file)
  }
  print(chi2_cytoRegions[order(p.corrected)])
  print(paste0("chi2_cytoRegions #rows: ", nrow(chi2_cytoRegions)))
  print(paste0("chi2_cytoRegions significant #rows: ", nrow(chi2_cytoRegions[chi2_cytoRegions$significant==T,])))
  
  #### BOXPLOTS & KW - Compute differences between pairs of groups
  if(nrow(chi2_cytoRegions!=0)){
    save_plot <- config.DiffMeth$overwrite_results
    gg <- plot_cytoRegions_lengthDistribution(cytoRegionsInfo, bins=30, config.DiffMeth, save_plot)
    # select only significant regions
    chi2_sig_cytoRegions <- chi2_cytoRegions[chi2_cytoRegions$significant==T,]
    #### HEATMAP of p.corrected values - only significant ones
    gg <- plot_heatmap_pCorrected(chi2_sig_cytoRegions, contexts=config.DiffMeth$contexts, region_type="CytoRegions", config.DiffMeth, save_plot)
    #### KW & DESCRIPTIVE STATS
    groupsCytoRegionsStats <- plot_boxplots_kw(methData_cytoRegions, cytoRegionsInfo, chi2_sig_cytoRegions, region_type="CytoRegions", config.DiffMeth, save_plot)
    groupsCytoRegionsStats$kw_results <- meth_data_add_main_cols(groupsCytoRegionsStats$kw_results, methData_cytoRegions)
    check_methData_cytoRegions(groupsCytoRegionsStats$kw_results)
    if(config.DiffMeth$overwrite_results){
      fwrite(groupsCytoRegionsStats$kw_results, file.path(config.DiffMeth$results_path, "groupsCytoRegionsKWStats.csv"))
      fwrite(groupsCytoRegionsStats$desc_stats, file.path(config.DiffMeth$results_path, "groupsCytoRegionsDescStats.csv"))
    }
  }
  
  ###################################################
  #### RNA EXPRESSION #### [Optional]
  if(config.DiffMeth$use_gene_expression_data){
    chi2_sig_cytoRegions <- chi2_cytoRegions[chi2_cytoRegions$significant==T,]
    #meanBeta_cytoRegions <- calc_mean_beta(methData_cytoRegions, chi2_cytoRegions[!is.na(chi2_cytoRegions$pval),], gene_expression_data)
    meanBeta_cytoRegions <- calc_mean_beta(methData_cytoRegions, chi2_cytoRegions, gene_expression_data)
    if(config.DiffMeth$overwrite_results){
      saveRDS(meanBeta_cytoRegions, file.path(config.DiffMeth$results_path, "meanBeta_cytoRegions.rds"))
    }
    # filter to the only significant regions in specific context
    meanBeta_cytoRegions <- meanBeta_cytoRegions[paste0(meanBeta_cytoRegions$context,"_",meanBeta_cytoRegions$region) %in% unique(paste0(chi2_sig_cytoRegions$context,"_",chi2_sig_cytoRegions$region))]
    #### PLOT HEATMAPS for all contexts and significant regions
    # cg <- split(meanBeta_cytoRegions, by="context", drop=TRUE)$CG
    # r <- chi2_sig_cytoRegions[chi2_sig_cytoRegions$context=="CG" & chi2_sig_cytoRegions$p.corrected<0.00001,"region"]
    # plot_region_sample_heatmap(cg[cg$region %in% r$region], "CytoRegions", config.DiffMeth, save_plot=F, plot_size=12, plot_margin=8, plot_gene=T)
    
    ret <- lapply(split(meanBeta_cytoRegions, by="context", drop=TRUE), 
                  function(x){plot_region_sample_heatmap(x, "CytoRegions", config.DiffMeth, save_plot, plot_size=24, plot_margin=8, plot_gene=T)})
    #### CALC slopeStatsRegions
    statsCytoRegions <- calcSlopeStats(meanBeta_cytoRegions, methData_cytoRegions)
    #select all sig_regions
    #selected_cytoRegion_context <- apply(unique(chi2_sig_cytoRegions[,c("region","context")]), 1, paste, collapse=".")
    selected_cytoRegion_context <- unique(apply(head(statsCytoRegions$slope[,c("region","context")],config.DiffMeth$top_beta_vs_exp_plots), 1, paste, collapse="."),
                                          apply(statsCytoRegions$corr[statsCytoRegions$corr$corr_pval<config.DiffMeth$max_pval,c("region","context")], 1, paste, collapse="."))
    #### PLOT BETA VS EXPR - general and detailed
    meanBeta_cytoRegions <- split(meanBeta_cytoRegions,list(meanBeta_cytoRegions$region, meanBeta_cytoRegions$context), drop=TRUE)
    ret <- lapply(meanBeta_cytoRegions[selected_cytoRegion_context], function(x){plot_beta_vs_exp(x, "CytoRegions", config.DiffMeth, save_plot)})
    statsCytoRegions <- statsRegions_rollup(statsCytoRegions, selected_cytoRegion_context)
    meanBeta_cytoRegions <- rbindlist(meanBeta_cytoRegions)
    if(config.DiffMeth$overwrite_results){
      fwrite(statsCytoRegions, file.path(config.DiffMeth$results_path, "slopeStatsCytoRegions.csv"))
    }
  }
    
}

