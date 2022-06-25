## ---------------------------

## Purpose of script: discover differentially methylated DNA regions
##
## Author: Michal Draminski
##
## Date Created: 2020-03-22
##

## ---------------------------


###################################################
#TODO to trzeba poprawic gdzie sortowac przycinac i filtrowac
#### KPI ####
if(config.DiffMeth$use_gene_expression_data){
  config.DiffMeth$text_size <- 8
  kpi_context <- "CG"
  config <- config.DiffMeth
  
  #### KPI REGIONS ####
  kpi_max_plot_regions <- 30
  kw_regions <- groupsRegionsStats$kw_results
  tests_number <- nrow(kw_regions[context==kpi_context,])
  kw_regions_plot <- kw_regions[significant==TRUE & context==kpi_context]
  kw_regions_plot <- head(kw_regions_plot, min(kpi_max_plot_regions, nrow(kw_regions_plot)))
  gg <- plot_regions_kpi(methData, meanBeta_regions, kw_regions_plot, kpi_context, region_main=NA, "Regions", tests_number, config.DiffMeth)
  file_name <- paste0(config$project_name,"_KPI_regions_",kpi_context,".",config$plot_format)
  ggsave(filename=file_name, plot = gg, device = config$plot_format, path = config$plots_path, dpi=120, width=3*nrow(kw_regions_plot), height=3*3, units = c("in"), limitsize = FALSE)
  
  #### KPI CYTOREGIONS ####
  kpi_max_plot_regions <- 10
  kw_cytoRegions <- groupsCytoRegionsStats$kw_results
  kpi_regions <- unique(kw_cytoRegions[significant==TRUE & context==kpi_context,]$region_main)
  #i=1
  for(i in 1:length(kpi_regions)){
    print(paste0("Region: ", kpi_regions[i]))
    tests_number <- nrow(kw_cytoRegions[context==kpi_context & region_main==kpi_regions[i],])
    kw_regions_plot <- kw_cytoRegions[significant==TRUE & context==kpi_context & region_main==kpi_regions[i]]
    kw_regions_plot <- head(kw_regions_plot, min(kpi_max_plot_regions, nrow(kw_regions_plot)))
    if(nrow(kw_regions_plot)>0){
      #order by position 
      setorder(kw_regions_plot,region)
      gg <- plot_regions_kpi(methData_cytoRegions, meanBeta_cytoRegions, kw_regions_plot, kpi_context, region_main=kpi_regions[i], "CytoRegions", tests_number, config.DiffMeth)
      file_name <- paste0(config$project_name,"_KPI_cytoRegions_",kpi_context,"_",kpi_regions[i],".",config$plot_format)
      ggsave(filename=file_name, plot = gg, device = config$plot_format, path = config$plots_path, dpi=120, width=3*nrow(kw_regions_plot), height=3*3, units = c("in"), limitsize = FALSE)
    }
  }  
}
