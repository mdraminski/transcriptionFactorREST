## ---------------------------

## Purpose of script: discover differentially methylated DNA regions
##
## Author: Michal Draminski
##
## Date Created: 2020-03-22
##

## ---------------------------

readDiffMethConfig <- function(file = "config.dm.yml"){
  cat(paste0("Reading the config file: '",file,"'\n"))
  ### READ MAIN CONFIG
  conf <- yaml::yaml.load_file(file)
  conf$project_path <- normalizePath(file.path(conf$project_path))
  if(conf$cyto_regions_find){
    conf$input_regions_path <- normalizePath(file.path(conf$input_regions_path))
  }
  conf$input_data_path <- normalizePath(file.path(conf$input_data_path))
  conf$input_data_def_path <- normalizePath(file.path(conf$input_data_def_path))
  if(conf$use_gene_expression_data){
    conf$input_expression_data_path <- normalizePath(file.path(conf$input_expression_data_path))
  }
  conf$meth_cut_poins <- c(-0.1, sort(conf$meth_cut_poins[conf$meth_cut_poins>0 & conf$meth_cut_poins<1]), 1)
  
  # Prepare local enviroment with project_enviroment()
  conf <- c(conf, prepareDiffMethEnviroment(project_path = conf$project_path, project_name = conf$project_name))
  conf$input_data_def <- load_inputDataDef(conf)

  return(conf)
}  

###################
# Function::prepareDiffMethEnviroment
###################
#it prepares the directories for analysis results
#project_path = "/home/mdraminski/TEMP2/"; project_name = "MGMTgene3";
prepareDiffMethEnviroment <- function(project_path, project_name){

  myconfig <- list(plots_path = file.path(project_path, project_name, "plots"),
                   results_path = file.path(project_path, project_name, "results"),
                   intersect_data_path = file.path(project_path, project_name, "intersect"),
                   tmp_path = file.path(project_path, project_name, "tmp"))
  
  tmp.res <- dir.create(file.path(project_path,project_name), showWarnings = F)
  lapply(myconfig, function(x){dir.create(x, showWarnings = F)})
  
  return(myconfig)
}

##################
#Function:: prepare_inputDataDef
##################
#input_data_path="/home/ZBO/FISHTANK/mdraminski/CytoMeth/results/methyl_results/"; save_def_file = F;
prepare_inputDataDef <- function(input_data_path, save_def_file = T){
  input_data_def <- data.frame(file_name = list.files(input_data_path, pattern = "methylation_results.bed"))
  input_data_def$sample_id <- sapply(strsplit(as.character(input_data_def$file_name), "[.]"), "[[", 1)
  input_data_def$group <- ""
  input_data_def$ignore <- "n"
  if(save_def_file){
    fwrite(input_data_def, file.path(input_data_path,"input_data_def.csv"))
  }
  return(input_data_def)
}

###################
#Function:: load_inputDataDef
##################
#config = config.DiffMeth
load_inputDataDef <- function(config){
  data_def <- as.data.frame(fread(config$input_data_def_path))
  check_inputDataDef(data_def,config)
  return(data_def)
}

###################
#Function:: check_inputDataDef
##################
#data_def = config.DiffMeth$input_data_def; config = config.DiffMeth
check_inputDataDef <- function(data_def, config){
  cat("Checking input_data_def file...\n")
  
  file_missing <- !file.exists(file.path(config$input_data_path, data_def$file_name))
  if(any(file_missing)){
    warning(paste0("Input Files are missing: ", paste0(file.path(config$input_data_path, data_def$file_name)[file_missing], collapse = "; ")))
  }
  expected_cols <- c("file_name","sample_id","group","ignore")
  cols_missing <- !expected_cols %in% names(data_def)
  if(any(cols_missing)){
    stop(paste0("Columns are missing: ", paste0(expected_cols[cols_missing], collapse = "; ")))
  }
  
  missing_group <- (is.na(data_def$group) | data_def$group=='') & data_def$ignore == "n"
  if(any(missing_group)){
    stop(paste0("Groups are missing for sample_id: ", paste0(data_def$sample_id[missing_group], collapse = "; ")))
  }
  
  data_def_groups <- table(data_def$group[!(is.na(data_def$group) | data_def$group=='') & data_def$ignore == "n"])
  if(length(data_def_groups)<2){
    stop(paste0("Two groups must be defined at least. Groups: ", paste0(names(data_def_groups), collapse = "; ")))
  }
  
  cat("File is correct.\n")
  return(TRUE)
}

###################
#Function:: add_TSS_promoter
##################
# start, end, strand columns must have names! Expected header of bed file: "chr","start","end","ENSG","HGNC","strand"

# this function depending on which strand is a gene selects the start column (+) or end (-) column as the TSS of a gene;
# if geneEnd = TRUE then another extra column is added; on + geneEnd = end; on - geneEnd = start
# if promoter=FALSE only TSS will be added if promoter=TRUE the promoter positions will be inserted
# the values of TSS positions will be changed to get the begining of the promoter and the end of the promoter;
# on + strand: start = TSS - upstream; end = TSS + downstream
# on - strand: start = TSS + upstream; end = TSS - downstream
# default upstream, downstream distancess == 2kb and could be changed
# because the start/end are replaced, new nolumnnames are inserted: prom_start / prom_end

#upstream = promoterUpstream; downstream = promoterDownstream; geneEnd=FALSE; promoter=TRUE;
add_TSS_promoter <- function(input_regions, upstream=2000, downstream=2000, geneEnd=FALSE, promoter=TRUE){
  expected_cols <- c("start","end","strand")
  if(any(!expected_cols %in% names(input_regions))){
    stop("Column names not as in the bed format: ",paste0(expected_cols, collapse = ','),"!")
  }

  input_regions$TSS[input_regions$strand=="+"] <- input_regions$start[input_regions$strand=="+"] #if + the TSS=start
  input_regions$TSS[input_regions$strand=="-"] <- input_regions$end[input_regions$strand=="-"] #if - the TSS=end
  
  if(geneEnd == TRUE){
    input_regions$geneEnd[input_regions$strand=="+"] <- input_regions$end[input_regions$strand=="+"] #if + the geneEnd=end
    input_regions$geneEnd[input_regions$strand=="-"] <- input_regions$start[input_regions$strand=="-"] #if - the geneEnd=start
  }
  
  if(promoter==TRUE){
    input_regions$start[input_regions$strand=="+"] <- input_regions$TSS[input_regions$strand=="+"] - upstream #if strand + decrease start
    input_regions$end[input_regions$strand=="+"] <- input_regions$TSS[input_regions$strand=="+"] + downstream #if strand + increase end
    input_regions$start[input_regions$strand=="-"]<-input_regions$TSS[input_regions$strand=="-"]- downstream  # this column keeps end of - strand gene promoter
    input_regions$end[input_regions$strand=="-"]<-input_regions$TSS[input_regions$strand=="-"] + upstream # this column keeps start of - strand gene promoter
  }
  input_regions$length <- input_regions$end - input_regions$start
  
  return(input_regions)
}

###################
# Function:: intersect_regions_with_samples
###################
# This function sorts the mask file and makes the intersection with CytoMethResults
#input_regions = inputRegions; config = config.DiffMeth;
intersect_input_data <- function(input_regions, config){
  input_regions <- setDT(input_regions)[,!names(inputRegions) %in% c("ensembl","hgnc"),, with=FALSE]
  output_folder <- config$intersect_data_path
  cat(paste0("Running bedtools intersect...\n"))
  
  if(config$intersect_strand_specific){
    bedtools_command = "bedtools intersect -wa -wb -s"
  }else{
    bedtools_command = "bedtools intersect -wa -wb"
  }
  
  # save input data to be able to use bedtools intersect
  input_regions_path <- file.path(config$tmp_path, "input_regions_tmp.bed")
  input_regions <- setDT(input_regions)[order(chr, start, end)]
  
  if(config$overwrite_results==T | !file.exists(input_regions_path)){
    fwrite(input_regions, file = input_regions_path, sep = "\t", scipen = getOption('scipen', 100), quote = F, col.names = F, row.names = F)
  }
  #in any casse we cn still use bedtools sorting
  #input_regions_path <- bedtools_sort(input_regions_path) # bedtools a file
  input_meth_data_path <- file.path(config$input_data_path, config$input_data_def$file_name[config$input_data_def$ignore=='n'])

  # intersect the a file with CytoMethFiles
  for(i in 1:length(input_meth_data_path)){
    curr_output_path <- file.path(output_folder, basename(input_meth_data_path[i]))
    bedtools_intersect(input_regions_path, input_meth_data_path[i], curr_output_path, bedtools_command, config$overwrite_results)
  }
}

###################
#Function:: bedtools_intersect
##################
#a_file_path = input_regions_path; b_file_path = input_meth_data_path; output_path = config$intersect_data_path; overwrite_results = FALSE;
bedtools_intersect <- function(a_file_path, b_file_path, output_path, bedtools_command = "bedtools intersect -wa -wb ", overwrite_results = FALSE){
  if(overwrite_results == FALSE & file.exists(output_path)){
    cat(paste0("File: ", output_path, " already exists! Skipping current intersection!\n"))
  }else{
    system_command <- paste0(bedtools_command, " -a ", a_file_path, " -b ", b_file_path, " > ", output_path)
    print(system_command)
    system(system_command)
  }
  print("Intersection is Done!")
}

###################
#Function:: bedSort
##################
# sort the file in bed format and return it in the same directory with string "_sorted.bed"
# file_directory - infile directory - file is sorted
# add_pattern - this will be added to the old file name, by default: "_sorted.bed"
bedtools_sort <- function(file_path, add_pattern="_sorted.bed"){
  out_dir <- paste0(gsub(".bed","",file_path),add_pattern)
  commandTobedtoolsSystem<-paste0("bedtools sort -i ",file_path, " > ", out_dir)
  system(commandTobedtoolsSystem)
  return(out_dir)
}

##################
#Function:: makeOneDataSet
##################
# Make a df from all intersected files and specify the coverage
#config = config.DiffMeth;
read_methData <- function(input_folder, config){
  cat(paste0("Loading data from: ", input_folder,"\n"))
  input_data_path <- list.files(path = input_folder, full.names = T)
  input_data <- lapply(input_data_path, fread)
  names(input_data) <- basename(input_data_path)
  input_data <- as.data.frame(rbindlist(input_data, idcol = "file_name"))
  #remove duplicated chr column
  if(all(input_data$V1 == input_data$V9)){
    input_data$V9 <- NULL
  }
  #move sample to the last position
  input_data <- input_data[,c(2:ncol(input_data),1)]
  input_data <- left_join(input_data, config$input_data_def[,c("file_name", "sample_id", "group")], by = "file_name")
  input_data$file_name <- NULL
  
  expected_col_names <- c("chr", "start", "end", "region", "gene", "strand", "TSS", "chr_c", "start_c", "end_c", "context", "beta", "strand_c", "cov", "numCs", "numTs", "position_c", "sample_id", "group")
  if(ncol(input_data) == 19){
    names(input_data) <- expected_col_names
  }else if(ncol(input_data) == 18){
    names(input_data) <- expected_col_names[expected_col_names!="TSS"]
  }else {
    print(paste0("The data contain ", ncol(input_data), " columns - 19 or 18 columns are expected!"))
    print(paste0("Expected columns [TSS is optional]: ", paste0(expected_col_names, collapse = ",")))
    print(head(input_data))
    return(input_data)
  }
  #make sure these are not factors
  input_data$sample_id <- as.character(input_data$sample_id)
  input_data$group <- as.character(input_data$group)
  
  input_data$chr_c <- NULL
  input_data <- input_data[input_data$cov >= config$min_coverage,]

  return(as.data.table(input_data))
}

##################
#Function:: read_inputRegions
read_inputRegions <- function(input_regions_path){
  input_regions <- fread(input_regions_path, header = FALSE)
  check_inputRegions(input_regions)
  names(input_regions) <- c("chr","start","end","region","gene","strand")
  input_regions$length <- input_regions$end - input_regions$start
  return(input_regions)
}

##################
#Function:: check_inputRegions
check_inputRegions <- function(input_regions){
  if(ncol(input_regions)!=6){
    warning(paste0("Error input_regions does not have 8 columns!"))
    print(input_regions)
    return(FALSE)
  }
  return(TRUE)
}

##################
#Function:: get_regions
##################
# meth_data=methData;
get_regions <- function(meth_data){
  #drop column sample_id to remove duplicated cytosines
  meth_data <- unique(meth_data[,c("chr","start","end","region","start_c","end_c","position_c")])
  regions <- meth_data[,.(cntC=.N), by=.(chr,start,end,region)]
  names(regions) <- c("chr","start","end","region","cntC")
  regions <- regions[,c("chr","start","end","region","cntC")]
  regions$length <- regions$end - regions$start
  regions$cyto_enrich <- regions$cntC/regions$length
  
  return(regions)
}










##################
#Function:: meth_data_remove_main_cols
##################
meth_data_remove_main_cols <- function(meth_data){
  meth_data <- unique(setDT(meth_data)[,!endsWith(names(meth_data),"_main"), with=FALSE])
  return(meth_data)
}

##################
#Function:: meth_data_add_main_cols
##################
#data_cytoRegions=groupsCytoRegionsStats$kw_results; meth_data_cytoRegions=methData_cytoRegions;
meth_data_add_main_cols <- function(data_cytoRegions, meth_data_cytoRegions){
  
  if(!any(c("chr","region") %in% names(data_cytoRegions))){  
    warning("CytoRegions data does not contain columns: 'chr' and 'region'")
    return(NULL)
  }
  data_cols <- names(data_cytoRegions)
  data_cytoRegions <- left_join(data_cytoRegions, unique(meth_data_cytoRegions[,c("chr","region","start_main","end_main","region_main")]), by=c("chr","region")) %>% data.table()
  region_idx <- which(data_cols=='region')
  new_cols <- c(data_cols[1:region_idx], names(data_cytoRegions)[endsWith(names(data_cytoRegions),"_main")], data_cols[(region_idx+1):length(data_cols)])
  data_cytoRegions <- setDT(data_cytoRegions)[,..new_cols]
  setorder(data_cytoRegions, kw_pval, na.last=T)
  
  return(data_cytoRegions)
}

##################
#Function:: get_groups_pairs
##################
get_groups_pairs <- function(meth_data){
  pairsToTest <- t(combn(unique(meth_data$group),2))
  pairsToTest <- setNames(data.frame(pairsToTest, stringsAsFactors = F), c("grp1","grp2"))
  pairsToTest$groupLabel <- paste0(pairsToTest$grp1, "_vs_", pairsToTest$grp2)
  return(pairsToTest)    
}

##################
#Function:: getContextDict
##################
getContextDict <- function(meth_data=NULL){
  context_dict <- list(allC = "", CG = c("CG"), nonCG = c("CHG", "CHH"))
  if(!is.null(meth_data)){
    context_dict[[which(context_dict=="")]] <- as.character(unique(meth_data$context))
  }
  return(context_dict)
}

##################
#Function:: meth_data_filter_context
##################
meth_data_filter_context <- function(meth_data, context=c("CG","nonCG","allC")){
  context <- context[1]
  context_values <- getContextDict(meth_data)
  context_values <- context_values[[context]]
  meth_data <- setDT(meth_data)[context %in% context_values,]
  return(meth_data)
}

##################
#Function:: meth_data_filter_regions
##################
meth_data_filter_regions <- function(meth_data, regions){
  meth_data <- setDT(meth_data)[region %in% unique(regions),]
  return(meth_data)
}
  
##################
#Function:: calcChi2MethDiff
##################
#meth_data=methData; contexts=c("allC","CG","nonCG"); config=config.DiffMeth;
#meth_data=methData_cytoRegions; contexts=c("allC","CG","nonCG"); config=config.DiffMeth;
calcChi2MethDiff <- function(meth_data, contexts=c("allC","CG","nonCG"), config){
  meth_data <- meth_data_remove_main_cols(meth_data)
  pairsToTest <- get_groups_pairs(meth_data)
  res_list <- list()
  #a loop over contexts i=1
  for(i in 1:length(contexts)){
    print(paste0("Context: ", contexts[i]))
    meth_data_context <- meth_data_filter_context(meth_data, contexts[i])
    #a loop over groups e.g. GIV vs GI j=1
    for(j in 1:nrow(pairsToTest)){
      groupsLabel <- pairsToTest[j,]$groupLabel #this is the string in the result list name
      print(paste0("Processing group label: ", groupsLabel, " [",j,"/",nrow(pairsToTest),"]"))
      meth_data_group <- meth_data_context[meth_data_context$group %in% c(pairsToTest[j,]$grp1, pairsToTest[j,]$grp2),]
      print(table(meth_data_group$group))
      current_ch2_result <- calcChi2StatsByRegions(meth_data_group, regionColumn = "region", geneColumn="gene", groupColumn="group", valueColumn = "beta", config)
      res_list_id <- paste0(groupsLabel, "_", contexts[i])
      if(nrow(current_ch2_result)>0){
        res_list[[res_list_id]] <- data.table(groups=groupsLabel, context=contexts[i], current_ch2_result, stringsAsFactors = F)
      }
    }
  }
  res_list <- rbindlist(res_list)

  return(res_list)
}

##################
#Function:: calcChi2StatsByRegions
##################
# FUNCTION DESCRIPTION
#input_data: dataset of methylations that result from intersection; it shall be a data frame
#regionColumn: this column vector unique values are used to run for loop in the run_test function: for a gene; a region; a promoter...
#filterColumn: this is a name of the columns where is the context of C: CpG, CpG, CHH...
#filterOutValue="" not samples from context column !CG, !CHH
#groupColumn: the name of the column where are the labels of the two compared groups
#twoLabelGroups: the strings of the two compared groups for which I run the test i.e. PA|HGG
#zeros: if Cs = 0 in all samples and present in PA and HGG - it was already removed from the function

#input_data=target_context_data; regionColumn="region"; geneColumn="gene"; groupColumn="group"; valueColumn = "beta";
calcChi2StatsByRegions <- function(input_data, regionColumn, geneColumn, groupColumn, valueColumn, config){
  cutPoints <- config$meth_cut_poins
  pval_max <- config$pval_max
  pval_corr <- config$pval_corr
  
  chi2Results <- list()
  input_data <- as.data.frame(input_data)
  region_gene_link <- unique(input_data[,c(regionColumn, geneColumn)])
  input_data <- split(setDT(input_data[,c(groupColumn,valueColumn,regionColumn)]), by=regionColumn)
  regions <- names(input_data)

  #run Chi2 test to select sigificantly different regions for is safer
  #chi2Results <- lapply(input_data, function(x) {data.frame(run_chi2_test(x, groupColumn, valueColumn, cutPoints))})
  for(i in 1:length(regions)){
      chi2Results[[i]] <- data.frame(run_chi2_test(input_data[[i]], groupColumn, valueColumn, cutPoints))
  }
  names(chi2Results) <- regions
  chi2Results <- rbindlist(chi2Results, idcol = regionColumn)
  
  if(nrow(chi2Results)>0){
    #chi2Results <- chi2Results[!is.na(chi2Results$pval) & !is.nan(chi2Results$pval) & chi2Results$pval < pval_max,]
    #correction for multiple testing
    chi2Results$p.corrected <- p.adjust(chi2Results$pval, method = pval_corr)
    #and add gene column
    chi2Results <- left_join(chi2Results, region_gene_link, by = regionColumn) %>% data.frame()
    chi2Results <- chi2Results[,c(regionColumn,geneColumn,"tstat","pval","p.corrected")]
    chi2Results$significant <- FALSE 
    chi2Results$significant[!is.na(chi2Results$pval) & !is.nan(chi2Results$pval) & chi2Results$p.corrected < pval_max] <- TRUE
    setDT(chi2Results)
  }
  return(chi2Results)
}

##################
#Function:: run_chi2_test
##################
#Function run_test performs Chi2 stat test on counts of beta values among defined intervals 
#(-0.25,0.2] (0.2,0.6] (0.6,1]
#GII/GIII        4223        98     169
#IDHmut          3552       144     213
#input_test_data = input_data[[1]];
run_chi2_test <- function(input_test_data, groupColumn, valueColumn, cutPoints){
  input_test_data <- as.data.frame(input_test_data)
  tbldata <- t(table(cut(input_test_data[,valueColumn], cutPoints), as.character(input_test_data[,groupColumn])))
  if(nrow(tbldata)==2){ #check if both groups are in the data
    test_res <- chisq.test(tbldata)
    test_res_df <- data.frame(tstat = test_res$statistic, pval = test_res$p.value)
    return(test_res_df)
  }else{
    return(NULL)
  }
}

##################
#Function:: heatMap_pCorrected
##################
# To make a heatmap showing for each region in C's context the p.corrected values
# chi2_results = chi2_results_regions; config = config.DiffMeth;
# chi2_results = chi2_results_cyto_regions; config = config.DiffMeth;
plot_heatmap_pCorrected <- function(chi2_results, contexts=c("allC","CG","nonCG"), region_type="Regions", config, save_plot=T){
  chi2_results <- chi2_results[context %in% contexts,]
  chi2_results$gr_con <- paste(chi2_results$context, chi2_results$groups, sep = "_")
  gg <- ggplot(chi2_results, aes(region, gr_con, fill=p.corrected)) + 
    geom_tile()+
    scale_fill_distiller()+
    theme_bw() +
    theme(text = element_text(size = config$text_size))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    labs(x = "Regions", y="Pairs in the C's context")

  if(save_plot){
    ggsave(filename=paste0(config$project_name, "_Heatmap_", region_type, ".",config$plot_format),
         plot = gg, device = config$plot_format, path = config$plots_path)
  }
  return(gg)
}

##################
#Function:: makeBoxplotsGroupsRegions
##################
# make plots for each significant region in context/contexts in which it was found significant
# meth_data = methData; regions_info=regionsInfo; region_type="Regions"; config = config.DiffMeth;
# meth_data = methData_kw; regions_info=cytoRegionsInfo; region_type="CytoRegions"; config = config.DiffMeth;
plot_boxplots_kw <- function(meth_data, regions_info, chi2_sig_regions, region_type="Regions", config, save_plot=T){
  if(nrow(chi2_sig_regions)==0){
    print("Significant regions data is empty!")
    return(NULL)
  }
  
  kw_results <- list()
  desc_stats <- list()
  contexts <- unique(chi2_sig_regions$context)
  #for over region_contexts i=1
  for(i in 1:length(contexts)){
    meth_data_kw <- meth_data_filter_context(meth_data, contexts[i])
    unique_regions <- unique(chi2_sig_regions$region[chi2_sig_regions$context==contexts[i]])
    #for loop - over regions j=1
    for(j in 1:length(unique_regions)){
      meth_data_region <- meth_data_kw[meth_data_kw$region==unique_regions[j],]
      kw_res <- plot_region_boxplots_kw(meth_data_region, context = contexts[i], region_type, tests_number=length(unique_regions), config)
      if(!is.null(kw_res)){
        if(kw_res$stats$significant & save_plot){
          ggsave(kw_res$file_name, plot = kw_res$gg, device = config$plot_format, path = config$plots_path)
        }
        desc_stats[[length(kw_results)+1]] <- kw_res$region_desc_stats
        kw_results[[length(kw_results)+1]] <- kw_res$stats
      }
    }
  }
  kw_results <- rbindlist(kw_results)

  #add more columns to the result
  final_cols <- c("chr","start","end","region","gene","context","cntC","length","cyto_enrich","kw_stat","kw_pval","kw_p.corrected","significant")
  kw_results <- left_join(kw_results, regions_info[,c("chr","start","end","region","cntC","length","cyto_enrich")], by=c("chr","region")) %>% data.table()
  kw_results <- kw_results[,..final_cols]
  setorder(kw_results, kw_pval, na.last=T)
  desc_stats <- rbindlist(desc_stats)
  desc_stats <- unique(desc_stats)
  setorder(desc_stats,chr,start,end,region,gene,context,group)
  resList <- list(kw_results=kw_results, desc_stats=desc_stats)
  return(resList)
}

##################
#Function:: plot_region_boxplots_kw
##################
#context = contexts[i]; tests_number=length(unique_regions);
plot_region_boxplots_kw <- function(meth_data_region, context, region_type="Regions", tests_number, config){
  meth_data_region <- meth_data_remove_main_cols(meth_data_region)
  meth_data_region <- meth_data_filter_context(meth_data_region, context)
  
  if(nrow(meth_data_region)==0){
    return(NULL)
  }
  if(length(unique(meth_data_region$region))>1){
    warning(paste0("One region is expected in meth_data_region. Found more!"))
    print(table(meth_data_region$region))
    return(NULL)
  }
  
  region <- unique(meth_data_region$region)[1]
  if(sum(table(meth_data_region$group)>0)>1){
    kw <- kruskal.test(beta ~ group, data = meth_data_region)
    region_desc_stats <- setDT(meth_data_region)[ ,.(min=min(beta,na.rm=T),q1=quantile(beta,na.rm=T)[2],median=median(beta,na.rm=T),mean=mean(beta,na.rm=T),q3=quantile(beta,na.rm=T)[3],max=max(beta,na.rm=T),n=.N),by=list(chr,start,end,region,gene,context,group)]
    meth_data_region <- as.data.frame(meth_data_region)
  }else{
    kw <- list(statistic=NaN, p.value=NaN)
    region_desc_stats = NULL
  }
  
  kw <- data.table(chr=meth_data_region$chr[1], start_region=meth_data_region$start_region[1], end_region=meth_data_region$end_region[1],
                   region=region, gene=meth_data_region$gene[1], context=context, kw_stat=kw$statistic, kw_pval=kw$p.value)
  kw$kw_p.corrected <- p.adjust(kw$kw_pval, method=config$pval_corr, n=tests_number)
  significant <- FALSE
  if(ifelse(is.nan(kw$kw_p.corrected),1,kw$kw_p.corrected) < config$pval_max){
    significant <- TRUE
  }
  kw <- cbind(kw, significant=significant)

  gg <- ggplot(meth_data_region, aes(x=group, y=beta, fill=group)) +
    geom_boxplot() +
    ggtitle(paste0(region," ",context, "\n(pval=",kw$kw_p.corrected,")")) +
    xlab("Subgroups of patients") +
    ylab("Methylation beta value") +
    theme_bw() +
    theme(text = element_text(size=config$text_size))+ 
    theme(plot.title = element_text(size=config$text_size)) +
    scale_fill_brewer(palette=config$brewer_palette) +
    coord_cartesian(ylim = c(0, 1)) +
    geom_jitter(shape=16, position=position_jitter(0.2))
    
  file_name <- paste0(config$project_name,"_Boxplot_",region_type,"_pval_",specify_decimal(kw$kw_p.corrected,config$precision),"_",unique(meth_data_region$region),"_", context,"_",config.DiffMeth$cyto_regions_max_distance, ".",config$plot_format)

  ret_val <- list(gg = gg, stats = kw, region_desc_stats = region_desc_stats, file_name = file_name)
  return(ret_val)
}

##################
#Function:: cyto_regions_length_distribution
##################
# To make a histogram plot of cyto_regions length
#cyto_regions=cytoRegions; config=config.DiffMeth;
#TODO MJD TO MA POPRAWIC ze bedzie piknie! mean/median i biny
plot_cytoRegions_lengthDistribution <- function(cyto_regions, bins=30, config, save_plot=T){
  gg <- ggplot(cyto_regions, aes(x=length)) +
    geom_histogram(aes(y=..density..), colour="black", fill="white", bins=bins)+
    geom_density(alpha=.2, fill="blue")+
    geom_vline(aes(xintercept = mean(length)),col='red',size=1)+
    geom_vline(aes(xintercept = median(length)),col='orange',size=1)+
    geom_text(x=15, y=0.07, label="mean", col="red", size = 6)+
    geom_text(x=15, y=0.06, label="median", col="orange", size=6)+
    theme_bw() +
    theme(text = element_text(size=config$text_size))+ 
    theme(plot.title = element_text(size=config$text_size+2))
    
  if(save_plot){
    ggsave(filename=paste0(config$project_name,"_Hist_CytoRegions_length_", config$cyto_regions_max_distance,".",config$plot_format),
         plot = gg, device = config$plot_format, path = config$plots_path)
  }
  return(gg)
}

##################
#RNA FUNCTIONS
##################

##################
#Function:: calc_mean_beta
##################
#meth_data = methData; chi2_sig_regions=chi2_regions[significant==T,]; gene_expression = gene_expression_data;
#meth_data = methData_cytoRegions; chi2_sig_regions = chi2_results_cyto_regions; gene_expression = gene_expression_data;
calc_mean_beta <- function(meth_data, chi2_sig_regions, gene_expression){
  meth_data <- meth_data_remove_main_cols(meth_data)
  mean_beta <- list()
  contexts <- unique(chi2_sig_regions$context)
  #for over region_contexts i=1
  for(i in 1:length(contexts)){
    cat(paste0("Context: ", contexts[i],"\n"))
    meth_data_filtered <- meth_data_filter_context(meth_data, contexts[i])
    sig_regions <- unique(chi2_sig_regions$region[chi2_sig_regions$context==contexts[i]])
    meth_data_filtered <- meth_data_filter_regions(meth_data_filtered, sig_regions)
    # for each region compute mean methylation
    mean_beta[[i]] <- meth_data_filtered[, .(mean_beta=mean(beta, na.rm=TRUE)), by=list(sample_id, region, gene, group)]
    mean_beta[[i]]$context <- contexts[i]
  }
  mean_beta <- rbindlist(mean_beta)
  mean_beta <- mean_beta[,c("context","region","gene","sample_id","group","mean_beta")]
  gene_expression <- setDT(gene_expression)[,c("sample_id", intersect(names(gene_expression), unique(mean_beta$gene))),with=FALSE]
  gene_expression <- melt(setDT(gene_expression),id.vars="sample_id",variable.factor = F)
  names(gene_expression) <- c("sample_id","gene","geneExp")
  setkey(mean_beta, sample_id, gene)
  setkey(gene_expression, sample_id, gene)
  mean_beta <- left_join(mean_beta, gene_expression, by=c("sample_id","gene"))
  
  return(mean_beta)
}

##################
#Function:: calcSlopeStats
##################
calcSlopeStats <- function(mean_beta_data, meth_data){
  mean_beta_data <- split(mean_beta_data,list(mean_beta_data$region, mean_beta_data$context), drop=TRUE)
  corr_stats <- rbindlist(lapply(mean_beta_data, calc_corr_params))
  setorder(corr_stats, corr_pval, na.last=T)
  slope_stats <- rbindlist(lapply(mean_beta_data, function(x){getModelParamsAll(x)}))
  slope_stats <- calc_slope_importance(slope_stats)
  meth_data_cols <- c("chr","start","end","region",names(meth_data)[endsWith(names(meth_data),"_main")],"gene")
  slope_stats <- right_join(unique(setDT(meth_data)[,..meth_data_cols]), slope_stats, by="region") %>% data.table()
  setorder(slope_stats, -Importance, na.last=T)
  retList <- list(corr=corr_stats,slope=slope_stats)
  return(retList)
}

##################
#Function:: plot_beta_vs_exp
##################
# For each region with computed mean meth
#mean_beta_context_region=meanBeta_regions[[1]]; config=config.DiffMeth;
plot_beta_vs_exp <- function(mean_beta_context_region, region_type="Regions", config, save_plot = T){
    geneName <- as.character(mean_beta_context_region$gene[1])
    current_context <- as.character(mean_beta_context_region$context[1])
    current_region <- as.character(mean_beta_context_region$region[1])
    importance <- calc_slope_importance(getModelParamsAll(mean_beta_context_region))
    importance <- round(importance$Importance[1],5)
    
    gg_detailed <- ggplot(mean_beta_context_region, aes(x=mean_beta, y=geneExp, color=group, shape=group))+
      geom_point(size=3)+
      #geom_smooth(method=lm, se=F, fullrange=T)+
      geom_smooth(method=lm, se=F)+
      xlab("Metylation beta value (mean)")+
      ylab(paste0(geneName," expression"))+
      #ggtitle(paste0("Region: ",current_region))+
      ggtitle(paste0(current_region))+
      theme_bw() +
      theme(text = element_text(size=config$text_size))+ 
      theme(plot.title = element_text(size=config$text_size)) +
      scale_color_brewer(palette=config$brewer_palette)
      
    gg_general <- ggplot(mean_beta_context_region, aes(x=mean_beta, y=geneExp))+
      geom_point(color="blue", shape=18, size=3)+
      #geom_smooth(method=lm, se=F, fullrange=T)+
      geom_smooth(method=lm, se=T)+
      xlab("Metylation beta value (mean)")+
      ylab(paste0(geneName," expression"))+
      #ggtitle(paste0("Region: ",current_region))+
      ggtitle(paste0(current_region," \n(corr=", specify_decimal(calc_corr_params(mean_beta_context_region)$corr_val,5),")"))+
      theme_bw() +
      theme(text = element_text(size=config$text_size))+ 
      theme(plot.title = element_text(size=config$text_size)) +
      scale_color_brewer(palette=config$brewer_palette)

    if(save_plot){
      ggsave(filename=paste0(config$project_name,"_BetaExpDetailed_",region_type,"_imp_",specify_decimal(importance,5),"_",current_region,"_",geneName,"_",current_context,".",config$plot_format),
             plot = gg_detailed, device = config$plot_format, path = config$plots_path)
  
      ggsave(filename=paste0(config$project_name,"_BetaExpGeneral_",region_type,"_imp_",specify_decimal(importance,5),"_",current_region,"_",geneName,"_",current_context, ".",config$plot_format),
             plot = gg_general, device = config$plot_format, path = config$plots_path)
    }
    
    retlist <- list(detailed=gg_detailed, general=gg_general)
    return(retlist)
}

##################
#Function:: plot_region_sample_heatmap
##################
#mean_beta=meanBeta_regions;
#mean_beta=split(meanBeta_cytoRegions, by="context", drop=TRUE)[[2]]; region_type="Regions"; config=config.DiffMeth;
plot_region_sample_heatmap <- function(mean_beta, region_type="Regions", config, save_plot=T, plot_size=12, plot_margin=6, plot_gene=F){
  if(nrow(mean_beta)==0)
    return(FALSE)
  
  current_context <- mean_beta$context[1]
  sample_label <- unique(mean_beta[,c("sample_id","group")])
  sample_label$group <- paste0("(",sample_label$group,")")
  sample_label$label <- apply(sample_label,1,function(x) paste(x, collapse = " "))
  heat_matrix <- as.data.frame(dcast(mean_beta, sample_id ~ region, value.var = "mean_beta"))
  rownames(heat_matrix) <- left_join(data.frame(sample_id=heat_matrix$sample_id),sample_label[,c("sample_id","label")], by="sample_id")$label
  heat_matrix$sample_id <- NULL
  heat_matrix <- as.matrix(heat_matrix)
  if(plot_gene){
    labels <- aggregate(gene ~ region, unique(mean_beta[,c("gene","region")]), paste)
    labels$label <- paste0(labels$region," (",labels$gene,")")
    colnames(heat_matrix) <- left_join(data.frame(region=colnames(heat_matrix)), labels, by="region")$label
  }
  if(region_type=="Regions"){
    margins<-c(plot_margin*2, 5)
  }else{
    margins<-c(plot_margin, 5)
  }
  file_name <- file.path(config$plots_path,paste0(config$project_name,"_Heatmap_beta_",region_type,"_",current_context, ".", config$plot_format))
  if(save_plot)
    openPlotFile(file_name, width = plot_size, height = plot_size, res = 72)
  heatmap(heat_matrix, scale="none", main=paste0("Context: ",current_context), xlab=region_type, margins=margins) #scale = c("row", "column")
  if(save_plot)
    dev.off()
  
  return(TRUE)
}

##############
#Function:: getModelParamsAll
#############
#data=mean_beta_context_region
getModelParamsAll <- function(data){
  data$group <- as.character(data$group)
  data_split_group <- split(data, data$group)
  fit_results <- lapply(data_split_group, function(x){lm(mean_beta ~ geneExp, data = x)})
  model_params <- lapply(fit_results, function(x) {getModelParams(x)})
  model_params <- rbindlist(model_params)
  model_params <- data.table(context=unique(data$context), region=unique(data$region), group=names(data_split_group), model_params)
  return(model_params)
}

#####################
#Function:: getModelParams
#####################
#lm_model=fit_results[[1]]
getModelParams <- function(lm_model){
  if(nrow(summary(lm_model)$coef)>1){
    pval <- round(summary(lm_model)$coef[2,4],5)
  }else{
    pval <- NA
  }
  df_lm_res <- data.table(AdjR2 = signif(summary(lm_model)$r.squared,5),
                          Intercept = signif(lm_model$coef[[1]],5),
                          Slope = signif(lm_model$coef[[2]],5),
                          pvalue = pval)
  return(df_lm_res)
}

#####################
#Function:: calc_slope_importance
#####################
calc_slope_importance <- function(slope_stats){
  slope_stats$Importance <- (slope_stats$Slope + abs(min(slope_stats$Slope, na.rm = T)))
  #boxplot(slopeStats$Importance)
  slope_stats_importance <- setDT(slope_stats)[,.(Importance=sd(Importance,na.rm=T)), by=list(context,region)][order(-Importance)]
  slope_stats_importance <- slope_stats_importance[,c("region","context","Importance")]
  return(slope_stats_importance)
}

#####################
#Function:: calc_corr_params
#####################
calc_corr_params <-function(x){
  if(nrow(x)<3){
    return (data.table())
  }else{
    corr_result <- cor.test(x$mean_beta,x$geneExp)
    corr_result <- data.table(context=x$context[1],region=x$region[1],gene=x$gene[1],corr_pval=corr_result$p.value,corr_val=corr_result$estimate)
    return(corr_result)
  }
}

#####################
#Function:: statsRegions_rollup
#####################
statsRegions_rollup <- function(stats_regions, selected_items){
  stats_regions <- left_join(stats_regions$slope, stats_regions$corr, by=c("context","region","gene"))
  stats_regions$saved <- FALSE
  stats_regions$saved[paste0(stats_regions$region,".",stats_regions$context) %in% selected_items] <- TRUE
  return(stats_regions)
}

##################
#Function:: plot_regions_kpi
##################
#meth_data=methData;mean_beta_data=meanBeta_regions;kw_results=groupsRegionsStats$kw_results;context="CG";region_main=NA;region_type="Regions";tests_number=NA;regions_limit=10;config=config.DiffMeth;
#meth_data=methData_cytoRegions;mean_beta_data=meanBeta_cytoRegions;kw_results=groupsCytoRegionsStats$kw_results;context="CG";config=config.DiffMeth;
#region_main <- kw_results[context=="CG" & significant==TRUE,]$region_main[1]
plot_regions_kpi <- function(meth_data, mean_beta_data, kw_results, context="CG", region_main=NA, region_type="Regions", tests_number=NA, config){
  curr_context <- context
  curr_region_main <- region_main
  meth_data <- meth_data_filter_context(meth_data, context=curr_context)
  mean_beta_data <- setDT(mean_beta_data)[context==curr_context,]

  if(is.na(curr_region_main)){
    if(is.na(tests_number))
      tests <- length(unique(kw_results$region[kw_results$context==curr_context]))
    sig_regions <- setDT(kw_results)[context==curr_context,]
  }else{
    if(is.na(tests_number))
      tests_number <- length(unique(kw_results$region[kw_results$context==curr_context & region_main %in% curr_region_main]))
    sig_regions <- setDT(kw_results)[context==curr_context & region_main %in% curr_region_main,]
  }
  
  regions_to_plot <- length(unique(sig_regions$region))
  print(paste0("Number of ",region_type," to plot: ",regions_to_plot))
  gg <- NULL
  if(regions_to_plot>0){
    gg1 <- list()
    gg2 <- list()
    gg3 <- list()
    for(i in 1:nrow(sig_regions)){
      meth_data_region <- meth_data[region==sig_regions$region[i],]
      gg <- plot_region_boxplots_kw(meth_data_region,curr_context,region_type=region_type,tests_number=tests_number,config)$gg
      if(i<nrow(sig_regions))
        gg <- gg + theme(legend.position = "none")
      if(i>1)
        gg <- gg + theme(axis.title.y=element_blank())
      gg1[[i]] <- gg
      tmpplots <- plot_beta_vs_exp(mean_beta_data[mean_beta_data$region==sig_regions$region[i],],region_type,config,save_plot=F)
      if(i<nrow(sig_regions))
        tmpplots$detailed <- tmpplots$detailed + theme(legend.position = "none")
      gg2[[i]] <- tmpplots$general
      gg3[[i]] <- tmpplots$detailed
    }
    gg <- wrap_plots(c(gg1,gg2,gg3),ncol=nrow(sig_regions),nrow=3)
  }  
  
  return(gg)
}

##################
#Function:: add_hgnc
##################
add_hgnc <- function(input_regions){
  require(GeneSummary)
  geneData <- loadGeneData()
  geneData <- geneData[,c("ensembl_gene_id","hgnc_symbol")]
  geneData <- unique(geneData)
  names(geneData) <- c("ensembl","hgnc")
  
  if(all(startsWith(toupper(input_regions$gene),"ENSG"))){
    #add hgnc
    input_regions$ensembl <- input_regions$region
    input_regions <- left_join(input_regions, geneData, by="ensembl")
  }else{
    #add ensg
    input_regions$hgnc <- input_regions$region
    input_regions <- left_join(input_regions, geneData, by="hgnc")
  }
  return(input_regions)
}

##################
#Function:: loadGeneData
##################
loadGeneData <- function(){
    ensembl <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
    # biomaRt::listAttributes(ensembl)
    # "phenotype_description"
    gene_data <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                               "hgnc_symbol",
                                               "chromosome_name",
                                               "start_position",
                                               "end_position",
                                               "description",
                                               "gene_biotype",
                                               "strand",
                                               "band",
                                               "transcript_count",
                                               "percentage_gene_gc_content"), mart = ensembl)
  return(gene_data)
}
