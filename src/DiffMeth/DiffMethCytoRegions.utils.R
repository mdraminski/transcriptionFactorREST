## ---------------------------

## Purpose of script: discover differentially methylated DNA regions
##
## Author: Michal Draminski
##
## Date Created: 2020-03-22
##

## ---------------------------

##################
#Function:: get_cytoRegions
##################
get_cytoRegions <- function(meth_data_cytoRegions){
  #recalc number of C in cyto_regions_info - keep in mind that position_c for strand + and - are bit different
  meth_data_cytoRegions$posC <- meth_data_cytoRegions$position_c
  meth_data_cytoRegions$posC[meth_data_cytoRegions$strand_c=="-"] <- meth_data_cytoRegions$posC[meth_data_cytoRegions$strand_c=="-"]-1
  cyto_regions_info <- meth_data_cytoRegions[,.(cntC=length(unique(posC))),by=list(chr,start,end,region)]
  meth_data_cytoRegions$posC <- NULL
  cyto_regions_info$length <- cyto_regions_info$end - cyto_regions_info$start
  cyto_regions_info$cyto_enrich <- cyto_regions_info$cntC/cyto_regions_info$length
  
  return(cyto_regions_info)
}

##################
#Function:: find_cytoRegions
##################
# inputData: the data based on which regions shall be estimated
# max_distance: what is the maximal distance between C-C
# meth_data=methData; max_distance=3;
find_cytoRegions <- function(meth_data, max_distance=10){
  meth_data <- split(meth_data, meth_data$chr)
  #using function:getCytoGroups assign each C to a group
  
  cat("Discovering C-rich regions...\n")
  #x <-meth_data[[1]]
  #unique(x[x$start_c>=153549899 & x$start_c<=153549912,c("chr","start_c","end_c","position_c")])
  #x[x$start_c==153549905,]
  
  cyto_regions <- lapply(meth_data, function(x) {
    x <- unique(x[,c("chr","start_c","end_c")])
    names(x) <-  c("chr", "start" , "end")
    x <- setDT(x)[order(chr, start, end)]
    res <- get_cyto_groups(x, max_distance)
    setnames(res, "group", "region")
    return(res)
  })
  
  #tu robimy zwijanie po grupach w ramach chromosomow
  cyto_regions <- lapply(cyto_regions, function(x){
    x <- setDT(x)[,.(start=min(start), end=max(end), cntC=.N, chr=unique(chr)), by=region]
    setDT(x)[order(region)]
  })
  
  #wrzucamy wszystko do jednego data frame
  cyto_regions <- rbindlist(cyto_regions)
  #no i jeszcze sie przyda kolumna ktora unikalnie rozroznia grupy bo to samo region moze byc w kilku chromosomach
  cyto_regions$region_label <- paste0(cyto_regions$chr,'_',cyto_regions$region)
  #to jeszcze kolumna z dlugoscia regionow
  cyto_regions$length <- cyto_regions$end - cyto_regions$start
  cyto_regions <- setDT(cyto_regions)[order(chr, start, end)]
  
  cyto_regions <- cyto_regions[,c("chr","start","end", "cntC", "region_label", "length")]
  setnames(cyto_regions, "region_label", "region")
  cyto_regions$cyto_enrich <- cyto_regions$cntC/cyto_regions$length
  
  return(cyto_regions)
}

##################
#Function:: get_Cyto_Groups
##################
get_cyto_groups <- function(x, distance){
  x <- as.data.table(x)
  x <- x[order(start)]
  x$diff <- c(1, diff(x$start, lag = 1))
  x$group <- x$diff
  x$group[x$group <= distance] <- 1
  x$group[x$group > distance] <- 0
  x$group_diff <- c(1,diff(x$group))
  x$group_diff[x$group_diff != 1] <- 0
  x$group_diff[x$group_diff == 1] <- 1:length(x$group_diff[x$group_diff == 1])
  x$group_diff_diff <- c(diff(x$group_diff, lag = 1),0)
  x$group_diff_diff[x$group_diff_diff<0] <- 0
  x$group_diff[x$group_diff_diff>0] <- x$group_diff_diff[x$group_diff_diff>0]
  x$group_diff_diff <- NULL
  x$group_diff[x$group_diff==0] <- NA
  x$group[!is.na(x$group_diff)]<-x$group_diff[!is.na(x$group_diff)]
  x$group_diff <- c(NA, na.omit(x$group_diff))[cumsum(!is.na(x$group_diff))+1]
  x$group[x$group == 1] <- x$group_diff[x$group == 1]
  x$group[x$group==0] <- NA
  x$group_diff <- NULL
  
  missing_group_num <- sum(is.na(x$group))
  if(missing_group_num>0){
    begin_id <- (max(x$group,na.rm = T)+1)
    x$group[is.na(x$group)] <- begin_id:(begin_id+missing_group_num-1)
    x$group <- as.numeric(factor(as.character(x$group), levels = unique(as.character(x$group))))
  }
  
  return(x)
}

##################
#Function:: cytoRegions_merge_methData
##################
#cyto_regions = cytoRegions; meth_data=methData; config = config.DiffMeth;
#cyto_regions = cytoRegions[cytoRegions$region=="chr17_36240",]; meth_data=methData; config = config.DiffMeth;
cytoRegions_merge_methData <- function(cyto_regions, meth_data, config){
  #Compute regions by chromosomes
  cyto_regions <- cyto_regions[,c("chr","start","end","region")]
  cyto_regions <- split(cyto_regions, cyto_regions$chr)
  
  if(config$intersect_strand_specific){
    meth_data <- meth_data[strand==strand_c,]
  }
  meth_data <- meth_data[,c("chr","start","end","region","gene","context","beta","strand_c","position_c","sample_id","group")]
  names(meth_data) <- c("chr","start_main","end_main","region_main","gene","context","beta","strand_c","position_c","sample_id","group")
  
  setkey(meth_data, "chr")
  #meth_data$region_main <- as.factor(meth_data$region_main)
  #meth_data$gene <- as.factor(meth_data$gene)
  #meth_data$context <- as.factor(meth_data$context)
  #meth_data$group <- as.factor(meth_data$group)
  meth_data <- split(meth_data, meth_data$chr)
  
  cat("Joining cyto regions and meth data. Plese wait...\n")
  for(i in 1:length(cyto_regions)){
    curr_chr <- names(cyto_regions)[i]
    cat(paste0("chr: ", names(cyto_regions)[i]))
    setkey(cyto_regions[[i]], "chr")
    if(!is.null(meth_data[[curr_chr]])){
      setkey(meth_data[[curr_chr]], "chr")
      cyto_regions[[i]] <- left_join_chr(cyto_regions[[i]], meth_data[[curr_chr]], bin_size = 1000)
    }else{
      cyto_regions[[i]] <- data.table()
    }
    cat("\n")
  }
  cat("Done!\n")
  cyto_regions <- rbindlist(cyto_regions) 
  return(cyto_regions)
}

##################
#Function:: left_join_chr
##################
#cyto_regions_dt=cyto_regions[[i]]; meth_data_dt=meth_data[[curr_chr]]; bin_size = 1000;
left_join_chr <- function(cyto_regions_dt, meth_data_dt, bin_size = 1000){
  merged_data <- list()
  setorder(cyto_regions_dt, "end")
  setkey(meth_data_dt, "position_c")
  if(bin_size>=nrow(cyto_regions_dt)){
    bin_idx <- c(1, nrow(cyto_regions_dt))
  }else{
    bin_idx <- c(0, seq(bin_size, nrow(cyto_regions_dt), bin_size))
    if(bin_idx[length(bin_idx)]<nrow(cyto_regions_dt))
      bin_idx <- c(bin_idx,nrow(cyto_regions_dt))
  }
  bins <- (length(bin_idx)-1)
  cat(paste0(" #bins: ", bins, " "))
  for(i in 1:bins){
    if(i==1){
      curr_idx <- (bin_idx[i]):bin_idx[i+1]
    }else{
      curr_idx <- (bin_idx[i]+1):bin_idx[i+1]
    }
    #cat(paste0("bin=",i," # ", min(curr_idx),":", max(curr_idx),"\n"))
    current_min_pos <- min(cyto_regions_dt$start[curr_idx])
    current_max_pos <- max(cyto_regions_dt$end[curr_idx])
    merged_data[[i]] <- left_join(cyto_regions_dt[curr_idx,], meth_data_dt[meth_data_dt$position_c>=current_min_pos & meth_data_dt$position_c<=current_max_pos,], by = "chr")
    merged_data[[i]] <- merged_data[[i]][start<=position_c & position_c<=end]
    merged_data[[i]] <- merged_data[[i]][start_main<=start & end<=end_main]
    cat(".")
  }
  merged_data <- rbindlist(merged_data)
  return(merged_data)
}

##################
#Function:: check_methData_cytoRegions
##################
check_methData_cytoRegions <- function(data_cytoRegions){
  broken_rows <- data_cytoRegions[!(data_cytoRegions$start>=data_cytoRegions$start_main & data_cytoRegions$end<=data_cytoRegions$end_main),]
  if(nrow(broken_rows)>0){
    warning("data_cytoRegions contains broken rows:")
    cat("Warning! data_cytoRegions contains broken rows:\n")
    print(broken_rows)
    return(FALSE)
  }
  return(TRUE)
}

