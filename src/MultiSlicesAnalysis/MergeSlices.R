##########----- 1. Load packages --##########
suppressPackageStartupMessages({
  library(argparse)
  library(stringr)
  library(dplyr)
  library(Seurat)
  library(SingleCellExperiment)
  library(ggplot2)
  library(patchwork)
  library(scater)
  library(harmony)
  library(BayesSpace)
  library(tictoc)
})
set.seed(149)


###########----- 2. Functions -----###########
##############################################
#     Main function 1: Load raw rds data     #
##############################################
#' @title Load rds data
#' @description Load raw rds data, transform Seurat object into SingleCellExperimrnt object, and combine all samples.
#' @name LoadData
#' @usage LoadData(chip_file, out, bin)
#'
#' @param chip_file A character vector specifying the names of files(comma-separated CSV format) containing chip information, each file contains 2 columns 'chip_id' and 'path'
#' @param out A string specifying the directory stored output files
#' @param bin An integer indicating the bin size
#'
#' @return A SingleCellExperiment object combined all sample information
#' @export
#'
LoadData <- function(chip_file, out, bin){
  options(warn = -1)
  ### Check parameters
  if(missing(out)) stop("Parameter 'out' should be provided!")
  if(! dir.exists(out)){
    dir.create(out, recursive = T)
    print(paste0("Output directory does not exist! Now create a output directory named '", out, "'"))
  }
  if(! is.numeric(bin)) stop("'bin' is not a numeric!")
  if(missing(chip_file)){
    stop("Parameter 'chip_file' should be provided!")
  }else{
    # Read chip file
    chipfile = read.csv(chip_file, header = T, sep = ",", stringsAsFactors = F)
  }
  
  
  ### Load rds data (Seurat --> SingleCellExperiment)
  # Set a list to store all samples
  tic("Load data")
  obj.list = list()
  for(i in 1:nrow(chipfile)){
    # check
    # print(chipfile$path[i])
    
    # Read rds
    seuratObj = readRDS(chipfile$path[i])
    
    # Get counts
    counts = as(seuratObj@assays$Spatial@data, "dgCMatrix")
    
    # Get colData
    df = data.frame(X = seuratObj@meta.data$x, Y = seuratObj@meta.data$y, stringsAsFactors = F)
    df$spot_id = paste(df$X, df$Y, sep = "X")
    colData = df[,c("spot_id","X","Y")]
    colnames(colData) = c("spot","row","col")
    rownames(colData) = colnames(counts)
    colData$sample_name = chipfile$chip_id[i]
    
    # Get rowData
    genes = rownames(counts)
    rowData = data.frame(genes, stringsAsFactors = F)
    
    # Create SingleCellExperiment object
    sce = SingleCellExperiment(assays = list(counts = counts),
                               rowData = rowData,
                               colData = colData)
    
    # Add in the list
    obj.list[chipfile$chip_id[i]] = sce
  }
  print(paste0("Load ", length(obj.list), " samples in a list!"))
  toc()
  
  
  ### Find features intersections of all samples to combine data
  tic("Find intersectional features")
  features = FindIntersetFeatures(obj.list = obj.list)
  toc()
  
  
  ### Combine all samples
  tic("Combine all samples")
  sce.combined = CombinedObj(obj.list = obj.list, features = features)
  print(paste0("Combine ", length(obj.list), " samples based on ", length(features), " intersectional features"))
  toc()
  
  
  ### Export metadata to transform coordinate
  meta = data.frame("spot" = sce.combined$spot,
                    "row" = sce.combined$row,
                    "col" = sce.combined$col,
                    "group" = sce.combined$sample_name, stringsAsFactors = F)
  rownames(meta) = colnames(sce.combined)
  write.csv(meta, paste0(out, "/sce_combined_meta.csv"), quote = F)
  
  
  ### Plot
  pdf(file = paste0(out, "/sce_combined_scatter_before_transcoord.pdf"), width = 10, height = 10)
  p <- plot(x = meta$row/bin, y = meta$col/bin, xlab = "", ylab = "")
  print(p)
  dev.off()
  
  
  saveRDS(object = sce.combined, file = paste0(out, "/sce_combined.rds"))
}

##############################################
#     Main function 1: Load raw rds data     #
#     Sub 1: Find features intersection      #
##############################################
#' @title Find features intersection of all samples
#' @description  Find features intersection of all samples to combine
#'
#' @param obj.list A SingleCellExperiment object list containing all sample information
#'
#' @name FindIntersetFeatures
#' @usage FindIntersetFeatures(obj.list)
#' @return A character vector containing features intersection of all samples
#' @export
#'
FindIntersetFeatures <- function(obj.list){
  # Set a list to store results
  gene.list = list()
  
  # Get all features of each sample in obj.list
  for(obj in 1:length(obj.list)){
    gene.list[[obj]] = rownames(obj.list[obj][[1]])
  }
  
  # Find intersection of all sample gene set
  features = Reduce(intersect, gene.list)
  
  # Print how many features found
  print(paste0(length(features), " intersectional features were found in ", length(obj.list), " samples"))
  
  return(features)
}

##############################################
#     Main function 1: Load raw rds data     #
#           Sub 2: Combine samples           #
##############################################
#' @title Combine samples
#' @description Combine all samples
#' @name CombinedObj
#' @usage CombinedObj(obj.list, features)
#'
#' @param obj.list A SingleCellExperiment object list containing all sample information
#' @param features A character vector containing features intersection of all samples
#'
#' @return A singleCellExperiment object after combine all samples
#' @export
#'
CombinedObj <- function(obj.list, features){
  # Filter each object based on features
  for(obj in names(obj.list)){
    obj.list[obj] = obj.list[obj][[1]][features,]
  }
  
  # Combine all filtered objects
  sce.combined = cbind(obj.list[[1]], obj.list[[2]], deparse.level = 2)
  if(length(obj.list) > 2){
    for(n in seq(3, length(obj.list))){
      sce.combined = cbind(sce.combined, obj.list[[n]], deparse.level = 2)
    }
  }
  
  return(sce.combined)
}


##############################################
# Main function 2: Integrate multiple slices #
##############################################
#' @title Integragte multiple slices
#' @description Integragte multiple slices, including coordinate transform, remove batch effects using harmony, cluster using bayesSpace
#' @name IntegrateMultiSlice
#' @usage IntegrateMultiSlice(sce, out, bin, transform, x_init_raw, y_init_raw, right, left, vertical, horizon, complete, num, n.PCs, q, nrep)
#'
#' @param sce The path of singleCellExperiment object after combine all samples
#' @param out A string specifying the directory stored output files
#' @param bin An integer indicating the bin size
#' @param transform Logical. Whether to transform coordinate? If TRUE, the default, transform coordinate
#' @param x_init_raw An integer indicating the starting x-coordinate
#' @param y_init_raw An integer indicating the starting y-coordinate
#' @param right A string specifying sample names need to turn right 90
#' @param left A string specifying sample names need to turn left 90
#' @param vertical A string specifying sample names need to reverse top and bottom
#' @param horizon A string specifying sample names need to reverse left and right
#' @param complete A string specifying sample names need to reverse top and bottom, left and right
#' @param num An integer indicating the sample number need to reset the starting x-coordinate
#' @param n.PCs An integer of principal components to compute, default: 50
#' @param q An integer of clusters
#' @param nrep An integer of MCMC iterations
#' @param pit.size Adjust point size for plotting
#' @param legend.size Adjust legend point size for plotting
#' @param legend.rows The desired number of rows of legends
#' @param width Set width of output figure
#' @param height Set height of output figure
#'
#' @return A singleCellExperiment object after clustering using bayesSpace
#' @export
#'
IntegrateMultiSlice <- function(sce, out, bin, transform = TRUE, x_init_raw = 0, y_init_raw = 0, right = NULL, left = NULL, vertical = NULL, horizon = NULL, complete = NULL, num = 8, n.PCs = 50, q = 16, nrep = 5000, pit.size = 1, legend.size = 3, legend.rows = 15, width = 10, height = 10){
  ### Check parameters
  if(missing(out)) stop("Parameter 'out' should be provided!")
  if(! dir.exists(out)){
    dir.create(out, recursive = T)
    print(paste0("Output directory does not exist! Now create a output directory named '", out, "'"))
  }
  if(! is.numeric(bin)) stop("'bin' is not a numeric!")
  if(! is.numeric(x_init_raw)) stop("'x_init_raw' is not a numeric!")
  if(! is.numeric(y_init_raw)) stop("'y_init_raw' is not a numeric!")
  if((! transform) & (length(c(right, left, vertical, horizon, complete)) > 0)) print("Flip coordinate only has an effect when 'transform' is TRUE!")
  if(! is.numeric(num)) stop("'num' is not a numeric!")
  if(! is.numeric(n.PCs)) stop("'n.PCs' is not a numeric!")
  if(! is.numeric(q)) stop("'q' is not a numeric!")
  if(! is.numeric(nrep)) stop("'nrep' is not a numeric!")
  if(! is.numeric(pit.size)) stop("'pit.size' is not a numeric!")
  if(! is.numeric(legend.size)) stop("'legend.size' is not a numeric!")
  if(! is.numeric(legend.row)) stop("'legend.row' is not a numeric!")
  if(! is.numeric(width)) stop("'width' is not a numeric!")
  if(! is.numeric(height)) stop("'height' is not a numeric!")
  
  
  ### Read combined object
  sce.combined = readRDS(sce)
  if(class(sce.combined) != "SingleCellExperiment") stop("'sce.combined' should be a SingleCellExperiment object!")
  
  
  #### Do coordinate transform if require
  if(transform){
    tic("Transform coordinate")
    print("Step1: Transform coordinate first!")
    # Get meta data to transform coordinate
    meta = data.frame("spot" = sce.combined$spot,
                      "row" = sce.combined$row/bin,
                      "col" = sce.combined$col/bin,
                      "group" = sce.combined$sample_name, stringsAsFactors = F)
    rownames(meta) = colnames(sce.combined)
    
    meta_adj = TransCoord(meta_data = meta,
                          x_init_raw = x_init_raw,
                          y_init_raw = y_init_raw,
                          right = right,
                          left = left,
                          vertical = vertical,
                          horizon = horizon,
                          complete = complete,
                          num = num)
    
    # Evaluate effect
    pdf(file = paste0(out, "/sce_combined_scatter_after_transcoord.pdf"), width = 10, height = 10)
    p1 <- plot(x = meta_adj$row_adj, y = meta_adj$col_adj, xlab = "", ylab = "")
    print(p1)
    dev.off()
    
    # Reset coordinate of combined object
    if(! identical(colnames(sce.combined), rownames(meta_adj))){
      meta_adj = meta_adj[match(colnames(sce.combined), rownames(meta_adj)),]
    }
    sce.combined@colData$spot = colnames(sce.combined)
    sce.combined@colData$row = meta_adj$row_adj
    sce.combined@colData$col = meta_adj$col_adj + abs(min(meta_adj$col_adj)) + 1
    sce.combined@colData$sample_name = meta_adj$group
    
    toc()
  }else{
    print("Step1: Don not need to transform coordinate!")
  }
  
  
  ### Batch correction
  print("Step2: Use Harmony to remove batch effect!")
  tic("Remove batch effect")
  # Preprocess
  sce.combined = spatialPreprocess(sce = sce.combined, n.PCs = n.PCs)
  
  # Check for batch effects first by plotting a low-dimensional representation of the data
  sce.combined = runUMAP(sce.combined, dimred = "PCA")
  colnames(reducedDim(sce.combined, "UMAP")) = c("UMAP1", "UMAP2")
  pdf(file = paste0(out, "/sce_combined_UMAP_before_removebatch.pdf"), width = 10, height = 10)
  
  p2 = ggplot(data = data.frame(reducedDim(sce.combined, "UMAP")),
              aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
    geom_point() +
    labs(color = "Sample") +
    theme_bw()
  
  print(p2)
  dev.off()
  
  # Use Harmony to integrate the samples
  sce.combined = RunHarmony(sce.combined, "sample_name", verbose = F)
  sce.combined = runUMAP(sce.combined, dimred = "HARMONY", name = "UMAP.HARMONY")
  colnames(reducedDim(sce.combined, "UMAP.HARMONY")) = c("UMAP1","UMAP2")
  
  pdf(file = paste0(out, "/sce_combined_UMAP_after_removebatch.pdf"), width = 10, height = 10)
  p3 = ggplot(data = data.frame(reducedDim(sce.combined, "UMAP.HARMONY")),
              aes(x = UMAP1, y = UMAP2, color = factor(sce.combined$sample_name))) +
    geom_point() +
    labs(color = "Sample") +
    theme_bw()
  
  print(p3)
  dev.off()
  
  toc()
  
  
  ### Use BayesSpace for clustering
  print("Step3: Use BayesSpace to cluster!")
  tic("Clustering")
  
  # cluster use BayesSpace
  sce.combined = spatialCluster(sce.combined, use.dimred = "HARMONY", q = q, nrep = nrep)
  saveRDS(object = sce.combined, file = paste0(out, "/sce_combined_bayes_cluster.rds"))
  
  # export coordinate data
  coord_df = as.data.frame(colData(sce.combined))
  write.csv(coord_df, file = paste0(out, "/sce_combined_bayes_cluster.csv"), quote = F, row.names = F)
  
  
  ### DimPlot
  colors = ColorPalette(number = length(unique(coord_df$spatial.cluster)))
  ClusterDimPlot(coord_df = coord_df, colors = colors, pit.size = 1, legend.size = 3, legend.rows = 15, width = 10, height = 10, out = out)
  
  
  ### Output a script of DimPlot
  script_file = paste0(out, "/DimPlot.R")
  
  # load packages
  cat("### Load packages\nlibrary(ggplot2)\nlibrary(RColorBrewer)\n\n", file = script_file, append = F)
  
  # functions
  cat("### Functions\n", file = script_file, append = T)
  cat("ColorPalette = ", file = script_file, append = T)
  cat(capture.output(dput(ColorPalette)), file = script_file, sep = "\n", append = T)
  cat("\nClusterDimPlot = ", file = script_file, append = T)
  cat(capture.output(dput(ClusterDimPlot)), file = script_file, sep = "\n", append = T)
  
  # load data and set colors
  cat("\n### Load data and set colors",
      "\ncoord_df = read.csv(file = paste0(out, '/sce_combined_adj_cluster_coordinate.csv'), header = T, stringsAsFactors = F)",
      "\ncolors = ColorPalette(length(unique(coord_df$spatial.cluster)))", file = script_file, append = T)
  
  # parameters
  cat("\n\n### Parameters",
      "\npit.size =", pit.size,
      "\nlegend.size =", legend.size,
      "\nlegend.rows =", legend.rows,
      "\nwidth =", width,
      "\nheight =", height,
      "\nout =", capture.output(dput(out)),
      file = script_file, append = T)
  
  # do DimPlot
  cat("\n\n### Do DimPlot",
      "ClusterDimPlot(coord_df = coord_df,
               colors = colors, pit.size = pit.size,
               legend.size = legend.size, legend.rows = legend.rows,
               width = width, height = height)",
      file = script_file, append = T, sep = "\n")
  
  toc()
}

##############################################
# Main function 2: Integrate multiple slices #
#        Sub 1: Transform coordinate         #
##############################################
#' @title Transform coordinate
#' @description Transform coordinate
#' @name TransCoord
#' @usage TransCoord(meta_data, x_init_raw = 0, y_init_raw = 0, right = NULL, left = NULL, vertical = NULL, horizon = NULL, complete = NULL, num = 8)
#'
#' @param meta_data A data.frame containing coordinate, 4 columns 'spot', 'row', 'col', 'group'
#' @param x_init_raw An integer indicating the starting x-coordinate
#' @param y_init_raw An integer indicating the starting x-coordinate
#' @param right A string specifying sample names need to turn right 90
#' @param left A string specifying sample names need to turn left 90
#' @param vertical A string specifying sample names need to reverse top and bottom
#' @param horizon  A string specifying sample names need to reverse left and right
#' @param complete A string specifying sample names need to reverse top and bottom, left and right
#' @param num An integer indicating the sample number need to reset the starting x-coordinate
#'
#' @return A data.frame containing transforming coordinate, 6 columns 'spot', 'row', 'col', 'group', 'row_adj', 'col_adj'
#' @export
#'
TransCoord <- function(meta_data, x_init_raw = 0, y_init_raw = 0, right = NULL, left = NULL, vertical = NULL, horizon = NULL, complete = NULL, num = 8){
  # Set a data.frame to store results
  results = data.frame()
  
  # Set start value
  x_init = 0
  y_init = 0
  x_offset = 50
  y_offset = 100
  x_max = 0
  x_min = 0
  y_max = 0
  y_min = 0
  x_dist = 0
  y_dist = 0
  
  # All sample names
  sample_name = unique(meta_data$group)
  
  
  # Flip transform list
  flip.list = list(right = right,
                   left = left,
                   vertical = vertical,
                   horizon = horizon,
                   complete = complete)
  
  # Do transform
  for(sample_n in 1:length(sample_name)){
    
    # Subset meta_data of sample_n
    sub_meta = subset(meta_data, group == sample_name[sample_n])
    
    # Check if the sample require flip transform
    idx = unlist(lapply(flip.list, FUN = function(x) sum(grepl(sample_name[sample_n], x))))
    types = names(idx[idx > 0])
    
    # Do flip transform if the sample require
    if(length(types) > 0){
      for(i in types){
        sub_meta = TransFlip(df = sub_meta, type = i)
      }
    }
    
    # Do coordinate transform
    x_min = min(sub_meta$row)
    y_max = max(sub_meta$col)
    
    gap_x = x_min - x_init
    gap_y = y_max - y_init
    
    sub_meta$row_adj = sub_meta$row - gap_x
    sub_meta$col_adj = sub_meta$col - gap_y
    
    x_max = max(sub_meta$row_adj)
    x_init = x_max + x_offset
    
    if(min(sub_meta$col_adj) < y_min){
      y_min = min(sub_meta$col_adj)
    }
    
    if(sample_n%%num == 0){
      x_init = x_init_raw
      y_init = y_min - x_offset
    }
    results = rbind(results, sub_meta)
  }
  return(results)
}

##############################################
# Main function 2: Integrate multiple slices #
#          Sub 2: Flip coordinate            #
##############################################
#' @title Flip coordinate
#' @description Flip coordinate, including turn right 90, turn left 90, reverse top and bottom, reverse left and right, reverse complete
#' @name TransFlip
#' @usage TransFlip(df, type)
#'
#' @param df A data.frame containing coordinate, 4 columns 'spot', 'row', 'col', 'group'
#' @param type 'right', 'left', 'vertical', 'horizon' or 'complete'. 'right': turn right 90; 'left': turn left 90; 'vertical': reverse top and bottom; 'horizon': reverse left and right; 'complete': reverse top and bottom, left and right
#'
#' @return A data.frame containing flip coordinate, 4 columns 'spot', 'row', 'col', 'group'
#' @export
#'
TransFlip <- function(df, type){
  # Check parameters
  types = c("right","left","vertical","horizon","complete")
  if(! type %in% types){
    stop("Type must be 'right', 'left', 'vertical', 'horizon' or 'complete' !")
  }else{
    tmp = df
    
    # 'right': Turn right 90, x = y, y = -x
    if(type == types[1]){
      tmp$row = df$col
      tmp$col = -df$row
      
      # 'left': Turn left 90, x = -y, y = x
    }else if(type == types[2]){
      tmp$row = -df$col
      tmp$col = df$row
      
      # 'vertical': Reverse top and bottom, x = x, y = -y
    }else if(type == types[3]){
      tmp$col = -df$col.
      
      # 'horizon': Reverse left and right, x = -x, y = y
    }else if(type == types[4]){
      tmp$row = -df$row
      
      # 'complete': Reverse top and bottom, left and right, x = -x, y = -y
    }else{
      tmp$row = -df$row
      tmp$col = -df$col
    }
  }
  
  return(tmp)
}

##############################################
# Main function 2: Integrate multiple slices #
#           Sub 3: Color palette             #
##############################################
#' @title  Color palette
#' @description Generate color palette
#' @name ColorPalette
#' @usage ColorPalette(number)
#'
#' @param number The number of colors used
#'
#' @return A string of colors
#' @export
#'
ColorPalette <- function(number){
  if(number < 25){
    colorScheme = c("#999999", "#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7", "#990000", "#9900cc", "#66FF66",
                    "#663300", "#0000FF", "#CC0033", "#FF0000", "#000099", "#660066",
                    "#333333", "#FFCCCC", "#993366", "#33CC33", "#000099", "#CC9900")
  }else{
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == "qual",]
    colorScheme = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
  }
  
  return(colorScheme)
}

##############################################
# Main function 2: Integrate multiple slices #
#              Sub 4: DimPlot                #
##############################################
#' @title DimPlot
#' @description DimPlot...
#' @name ClusterDimPlot
#' @usage ClusterDimPlot(coord_df, colors = NULL, pit.size = 1, legend.size = 3, legend.rows = 15, legend.title = "Clusters")
#'
#' @param coord_df A data frame of coordinate after cluster, must contain 3 columns: "row", "col" and "spatial.cluster"
#' @param colors A character of colors used
#' @param pit.size Adjust point size for plotting
#' @param legend.size Adjust legend point size for plotting
#' @param width Set width of output figure
#' @param height Set height of output figure
#' @param out A string specifying the directory stored output files
#' @param legend.rows The desired number of rows of legends
#'
#' @return A figure of DimPlot
#' @export
#'
ClusterDimPlot <- function(coord_df, colors = NULL, pit.size = 1, legend.size = 3, legend.rows = 15, width = 10, height = 10, out){
  pdf(file = paste0(out, "/sce_combined_DimPlot_bayes_cluster.pdf"), width = width, height = height)
  p = ggplot(data = coord_df, aes(x = row, y = col, color = as.factor(spatial.cluster))) +
    geom_point(shape = 19, size = pit.size) +
    coord_fixed() +
    labs(x = "", y = "") +
    guides(colour = guide_legend(override.aes = list(size = legend.size), nrow = legend.rows, title = "Clusters")) +
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "right") +
    theme_void()
  
  if(! is.null(colors)){
    p = p + scale_color_manual(values = colors)
  }
  
  print(p)
  dev.off()
}


###########----- 3. Do  work -----########
### do a test
# Load data
sce.combined = LoadData(chip_file = "/jdfssz1/ST_TSCBI/PROJECT_temp/USER/zhaodanhui/Notebook/test/input/chip_new.csv", out = "/jdfssz1/ST_TSCBI/PROJECT_temp/USER/zhaodanhui/Notebook/test/out_new2", bin = 50)

# Integrate multiple slices
right <- c('FP200000571BR_C5_2','FP200000571BR_A2_3', 'SS200000103BR_C4_2')
left <- c('FP200000571BR_A2_1','SS200000103BR_C4_1', 'FP200000571BR_C1_1',
          'FP200000571BR_C1_2', 'FP200000571BR_C5_1')
complete <- c('FP200000571BR_B1_2', 'FP200000571BR_B2_1', 'FP200000609BR_B5_1')

sce.combined = IntegrateMultiSlice(sce = "/jdfssz1/ST_TSCBI/PROJECT_temp/USER/zhaodanhui/Notebook/test/out_new2/sce_combined.rds", out = "/jdfssz1/ST_TSCBI/PROJECT_temp/USER/zhaodanhui/Notebook/test/out_new2", bin = 50, 
                                   transform = TRUE, x_init_raw = 0, y_init_raw = 0, 
                                   right = right, left = left, vertical = NULL, horizon = NULL, 
                                   complete = complete, num = 8, n.PCs = 50, q = 16, nrep = 5000, width = 10, height = 10)