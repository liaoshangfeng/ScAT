## Annotation analysis by SingleR

parser = argparse::ArgumentParser(description = 'Script for annotation analysis')
parser$add_argument('-i', '--input', dest = 'input', default=NULL, help = 'Input normalized data (*.txt.gz) [default %(default)s]')
parser$add_argument('-c', '--cluster', dest = 'cluster', default=NULL, help = 'Cluster type file name [default %(default)s]')
parser$add_argument('-o', '--out', dest = 'out', default=NULL, help = 'Directory to save file [default %(default)s]')
parser$add_argument('--prefix', dest = 'prefix', default = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident [default %(default)s]')
parser$add_argument('--ref',  dest = 'ref', default=NULL, help = 'Ref database name [default %(default)s]')
parser$add_argument('--ref_col',  dest = 'ref_col', default=NULL, help = 'Ref database col name [default %(default)s]')

opts = parser$parse_args()
print(opts)

## Debug


# opts <- list()
# opts$input <- "/home/tonyleao/wkd/data_to_viz/cell_annotation/data/e15_raw_matri_norm.txt.gz"
# opts$out <- "/home/tonyleao/wkd/data_to_viz/cell_annotation/out"
# opts$ref_col <- "/home/tonyleao/wkd/Scat/docs/celldex/HumanPrimaryCellAtlasData_coldata.rds"
# opts$ref <- "/home/tonyleao/wkd/Scat/docs/celldex/HumanPrimaryCellAtlasData.rds"
# opts$cluster <- "/home/tonyleao/wkd/data_to_viz/cell_annotation/data/e15_raw_matri_norm_clusterInfo.txt"
# print(opts)


### Function

#' ColorPalette
ColorPalette <- function(number){
  print(number)
  if (number <= 25){
    colorScheme = c("#999999", "#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7", "#990000", "#9900cc", "#66FF66",
                    "#663300", "#0000FF", "#CC0033", "#FF0000", "#000099", "#660066",
                    "#333333", "#FFCCCC", "#993366", "#33CC33", "#000099", "#CC9900")
  } else {
    suppressMessages(library(RColorBrewer))
    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    colorScheme <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
    return (colorScheme)
  }
}

#' define a function to remove the NA in SE object
.rm_NAs <- function(mat, rm.NA = "rows"){
    # Identify them first before removal, to ensure that
    # the same rows are columns are removed with 'both'.
    if(rm.NA == "rows" || rm.NA == "both"){
        keep_rows <- !rowAnyNAs(DelayedArray(mat))
    } else {
        keep_rows <- !logical(nrow(mat))
    }

    if (rm.NA=="cols" || rm.NA == "both") {
        keep_cols <- !colAnyNAs(DelayedArray(mat))
    } else {
        keep_cols <- !logical(ncol(mat))
    }
    # Avoid making unnecessary copies if possible.
    if (!all(keep_rows)) {
        mat <- mat[keep_rows,,drop=FALSE]
    }
    if (!all(keep_cols)) {
        mat <- mat[,keep_cols,drop=FALSE]
    }
    mat
}



#' define a function to perform cell annotation
#' tutorial: https://www.jianshu.com/p/9dcf1c0a4cfd
cell_annotation <- function(opts){
    data <- read.table(gzfile(opts$input), header = T, row.names = 1 , check.names = FALSE)
    
    # make SE object
    all.assays <- list()
    reference <- readRDS(opts$ref)
    all.assays[['logcounts']] <- .rm_NAs(reference, c("rows","cols","both","none"))

    arg <- list()
    ref_colname <- readRDS(opts$ref_col)
    arg$colData <- ref_colname
    # make the final SE object
    ref_se <- do.call(SummarizedExperiment, c(list(assays=all.assays),arg))

    pred1 <- SingleR(test = data, ref = ref_se, assay.type.test=1, labels = ref_se$label.main)
    pred2 <- SingleR(test = data, ref = ref_se, assay.type.test=1, labels = ref_se$label.fine)
    pred<- cbind(rownames(pred1), pred1$labels, pred2$labels)
    colnames(pred)<-c('cells','main','fine')
    write.table(pred, paste(opts$out, '/', opts$prefix, '_cell_types.txt', sep = ''), sep='\t', quote=FALSE, row.names = FALSE)
    return(pred)
}


#' define a function for cell component analysis
cell_component <- function(opts, pred){
    cell_type <- as.data.frame(pred)
    rownames(cell_type) <- cell_type[,1]
    cell_type <- cell_type[,-1]
    colnames(cell_type) = c("main", "fine")
    
    # 01 load cell cluster file
    cell_cluster <- read.table(opts$cluster, row.names = 1, header = T, check.names = F)
   
    # 02 merge the celltype and clusterInfo
    anno_cell <- intersect(rownames(cell_type), rownames(cell_cluster))
    df <- cell_type[rownames(cell_type) %in% anno_cell, ]
    
    df$cluster <- rep(NA, nrow(df))
    for (cell in anno_cell) {
      df[cell, ]$cluster  = cell_cluster[cell, ]
    }
  
    # 03 save the cluster cell component file
    write.table(df, paste(opts$out, '/', opts$prefix, '_cluter_cell_component.txt', sep = ''), sep = '\t', quote = FALSE, row.names = FALSE)

    #04  Draw a cell compnent histogram
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggthemes))
    df$number <- 1
    colorScheme <- ColorPalette(as.numeric(length(unique(df$main))))
    p <- ggplot(df, aes(cluster, number, fill = main)) +
      geom_bar(stat = "identity", position = "stack") +
      # coord_flip() +
      theme_bw() +
      # scale_fill_wsj("colors6", "")+
      scale_fill_manual(values = colorScheme[1:length(unique(df$main))]) +
      theme(axis.ticks.length = unit(0.5, 'cm')) +
      guides(fill = guide_legend(title = NULL))
    
    p <- p + theme(
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black")) + 
      labs(x = "Cell Number", y = "Cell Clsuter", title = "") + 
      theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(hjust = 1, vjust = -40, size = 14)
    )
    ggsave(paste(opts$out, '/', opts$prefix, '_cluster_cell_component.pdf', sep=''), plot = p, width = 0.5*length(unique(df$main)), height = 5)
    
  
    # 05 Draw a percentage stack diagram
    suppressMessages(library(plyr))
    df <- ddply(df, "cluster", transform,
                percent_weight = number / sum(number) * 100)
    
    p1 <- ggplot(df, aes(cluster, percent_weight, fill = main)) +
      geom_bar(stat = "identity", position = "stack") +
      # coord_flip() +
      theme_bw() +
      # scale_fill_wsj("colors6", "")+
      scale_fill_manual(values = colorScheme[1:length(unique(df$main))]) +
      theme(axis.ticks.length = unit(0.5, 'cm')) +
      guides(fill = guide_legend(title = NULL))
    
    
    p1 <- p1 + theme(
      axis.text.x = element_text(size = 14, color = "black"),
      axis.text.y = element_text(size = 14, color = "black")) + 
      labs(x = "Cell Cluster", y = "Cell component", title = "") + 
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(hjust = 1,vjust = -40,size = 14))
    ggsave(paste(opts$out, '/', opts$prefix, '_cluster_cell_component_percentage.pdf', sep=''), plot = p1, width = 0.5*length(unique(df$main)), height = 5)
}


#' define the function to perform workflow
main <- function(opts){
    knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
    suppressMessages(library(BiocStyle))
    suppressMessages(library(SingleR))
    suppressMessages(library(ExperimentHub))
    suppressMessages(library(DelayedArray))
    pred <- cell_annotation(opts)
    if (! is.null(opts$cluster)){
    cell_component(opts, pred)}
}




 ### Main Work flow
if ((! is.null(opts$input)) & (! is.null(opts$ref)) & (! is.null(opts$ref_col)) & (!is.null(opts$out))){
  if (! file.exists(opts$out)){
    dir.create(opts$out, recursive = TRUE)}
    main(opts)
  }else{
   parser$print_help()
   stop("You need to provide the following parameters: --input, --ref, --ref_col, --out", call.=FALSE)
}

