##########
#   Script for quality control and preprocessing analysis of stereo gem data. 
#   Function 1: load stereo-seq format matrix, convert to seurat object, and preview the basic qc metrics
#   Function 2: load seurat object to perform filteration and normalization
#   Note: 
#       filteration based on the maxFeature, minFeature, mitochondria UMI ratio per bin.
#       normalization based on SCTransfrom
#       If the default assays is 'Spatial', filtering and standardization will be performed by default, otherwise an error will be reported          
#           
#####
##2021/11/03
##Update: 
#   1.Support input *.gem.gz file
#   2.Add --fileType to distinguish *.matrix.txt.gz and *.gem.gz
##2022/02/07
##Update:
##  1.Support load *gem.gz file with different format, including 'MIDCounts', 'MIDCount', 'UMIDCount'
#####

### Get the parameters
library(argparse, quietly = TRUE)
parser = argparse::ArgumentParser(description = 'Script for quality control analysis of stereo-seq matrix')
parser$add_argument('-i', '--input', dest = 'input', help = 'input stereo-seq *.gem or seurat object *.rds filename')
parser$add_argument('-o', '--out', dest = 'out', help = 'directory where to save the output files, all output files will be indexed by sample ID')
parser$add_argument('-b', '--binsize', dest = 'binsize', default = 1, type = 'integer', help = 'bin size to binning, your input should be in binSize1 if you set this [default %(default)s]')

parser$add_argument('--ptsize', dest = 'ptsize', default = 1, type = 'double', help = 'point size for spatial dim plot [default %(default)s]')
parser$add_argument('--prefix', dest = 'prefix', default = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident [default %(default)s]')
parser$add_argument('--geneColumn', dest = 'geneColumn', default = 2, type = 'integer', help = 'Specify which column of genes.tsv or features.tsv to use for gene names [default %(default)s]')
parser$add_argument('--fileType', dest = 'fileType', default = 'gem', help = 'Set as matrix to support *.matrix.txt.gz, otherwise, no need to change  [default %(default)s]')


# basic qc metrics
parser$add_argument('--mtPattern', dest = 'mtPattern', default='^MT-', help = 'mitochondria gene pattern [default %(default)s]')
parser$add_argument('--rbPattern', dest = 'rbPattern', default='^RP[SL]', help = 'ribosome gene pattern [default %(default)s]')
parser$add_argument('--hbPattern', dest = 'hbPattern', default='^HB[^(P)]', help = 'hemoglobin gene pattern [default %(default)s]')
parser$add_argument('--platPattern', dest = 'platPattern', default='PECAM1|PF4', help = 'platelet gene pattern [default %(default)s]')

# filteration parms
# Currently, the minFeature, maxFeature and mtRatio were used for filteration
parser$add_argument('--mFilter', dest = 'mFilter', action="store_true", default=FALSE, help = 'Allow filteration analysis [default %(default)s]')
parser$add_argument('--minCount', dest = 'minCount', default = 0, type = 'integer', help = 'minimum UMI number [default %(default)s]')
parser$add_argument('--maxCount', dest = 'maxCount', default = 1000000, type = 'integer', help = 'maximum UMI number [default %(default)s]')
parser$add_argument('--minFeature', dest = 'minFeature', default = 0, type = 'integer', help = 'minimum Feature number [default %(default)s]')
parser$add_argument('--maxFeature', dest = 'maxFeature', default = 100000, type = 'integer', help = 'maximum Feature number [default %(default)s]')
parser$add_argument('--minCell', dest = 'minCell', default = 3, type = 'integer', help = 'minimum cell or BIN number express the gene (filter genes) [default %(default)s]')

parser$add_argument('--mtRatio', dest = 'mtRatio', default = 100, type = 'double', help = 'maximum mitochondria gene UMI ratio in a BIN [default %(default)s]')
parser$add_argument('--rbRatio', dest = 'rbRatio', default = 100, type = 'double', help = 'maximum ribosome gene UMI ratio in a BIN [default %(default)s]')
parser$add_argument('--hbRatio', dest = 'hbRatio', default = 100, type = 'double', help = 'maximum hemoglobin gene UMI ratio in a BIN [default %(default)s]')
parser$add_argument('--platRatio', dest = 'platRatio', default = 100, type = 'double', help = 'maximum platelet gene UMI ratio in a BIN [default %(default)s]')

opts = parser$parse_args()
print(opts)



## Functions

#' define function to transform *.gem file to seurat object
gem_to_seuratObject <- function(data, prefix = 'sample', binsize = 1){
    #' read the gem file  
    data <- fread(file = data)
    
    #' group counts into bins
    data$x <- trunc(data$x / opts$binsize) * opts$binsize
    data$y <- trunc(data$y / opts$binsize) * opts$binsize
    #trunc取整函数，只取整数部分
    #将横纵坐标根据binsize取整，比如，如果binszie=10,x=2401,y=4522,
    #那么处理后，x=2400,y=4520
    
    if ('MIDCounts' %in% colnames(data)) {
      data <- data[, .(counts=sum(MIDCounts)), by = .(geneID, x, y)]
    } else if ('UMICount' %in% colnames(data)) {
      data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
    } else if ('MIDCount' %in% colnames(data)) {
      data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]
    } 

    #' create sparse matrix from stereo
    data$cell <- paste0(prefix, ':', data$x, '-', data$y)
    data$geneIdx <- match(data$geneID, unique(data$geneID))
    data$cellIdx <- match(data$cell, unique(data$cell))
    
    mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                        dimnames = list(unique(data$geneID), unique(data$cell)))
    #sparseMatrix稀疏矩阵函数
    
    cell_coords <- unique(data[, c('cell', 'x', 'y')])
    #unique去除重复函数，删除cell,x和y都一样的行
    
    rownames(cell_coords) <- cell_coords$cell
    
    
    seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial', 
                                            names.delim = ':', meta.data = cell_coords)
    
    
    #' create pseudo image
    cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
    cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
    
    tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
    #matrix(aa,x,y)以aa为输入向量，创建一个x行y列的矩阵
    #构造一个seruat image
    
    tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                        tissue = 1,
                                        row = cell_coords$y, col = cell_coords$x,
                                        imagerow = cell_coords$y, imagecol = cell_coords$x)
    
    
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                     tissue_hires_scalef = 1,
                                     tissue_lowres_scalef = 1))
    #toJSON: 把json格式 转换成 list格式
    
    
    #' function to create image object
    generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
      if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
      }
      
      unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
      
      spot.radius <- unnormalized.radius / max(dim(x = image))
      
      return(new(Class = 'VisiumV1', 
                 image = image, 
                 scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                              fiducial = scale.factors$fiducial_diameter_fullres, 
                                              hires = scale.factors$tissue_hires_scalef, 
                                              lowres = scale.factors$tissue_lowres_scalef), 
                 coordinates = tissue.positions, 
                 spot.radius = spot.radius))
    }
    
    visiumV1Obj <- generate_spatialObj(image = tissue_lowres_image, 
                                       scale.factors = fromJSON(scalefactors_json), 
                                       tissue.positions = tissue_positions_list)
    #可以理解为构建一个spatial背景
    
    #' import image into seurat object
    visiumV1Obj <- visiumV1Obj[Cells(x = seurat_spatialObj)]
    DefaultAssay(visiumV1Obj) <- 'Spatial'
    
    seurat_spatialObj[['slice1']] <- visiumV1Obj
    rm(visiumV1Obj, data, mat)
    gc()
    return(seurat_spatialObj)
}


#' preview the data in Spatial
# 这个应该设置成数据预览的函数
preview_seuratObject <- function(sce, opts){
    # assay ('Spatial', 'RNA', 'SCT')
    nFeature <- paste('nFeature_', opts$assays, sep='')
    nCount <- paste('nCount_', opts$assays, sep='')

    # 图片的排版有问题，暂时分开输出
    # 01 plot QC metric with Violin plot
    pdf(file = paste(opts$out, '/', opts$prefix, '_QC_metrics_violinPlot.pdf', sep = ''), w = 10, h = 8)
    p <- suppressWarnings(VlnPlot(sce, group.by = 'orig.ident', features = c(nFeature, nCount, 'percent.mt', 'percent.ribo', 'percent.hb', 'percent.plat'), pt.size = 0.1, ncol = 3)) + NoLegend()
    print (p)
    # while (!is.null(dev.list()))  dev.off() 
    dev.off()
    
    # 02 plot featureScatter plot
    plot1 <- suppressWarnings(FeatureScatter(sce, feature1 = nCount, feature2 = 'percent.mt'))
    plot2 <- FeatureScatter(sce, feature1 = nCount, feature2 = nFeature)
    pdf(file = paste(opts$out, '/', opts$prefix, '_FeatureScatter.pdf', sep = ''), w = 8, h = 4)
    print (plot1 + plot2)
    # while (!is.null(dev.list()))  dev.off() 
    dev.off()

    # opts$ptsize = 10
    # if (F){
    if (opts$assays == 'Spatial'){
        # 03 plot nCount_Spatial plot
        pdf(file = paste(opts$out,'/', opts$prefix, '_nCount_Spatial',  '.pdf', sep = ''), w = 10, h = 5)
        p1  <- suppressWarnings(VlnPlot(sce, group.by = 'orig.ident', features = nCount, pt.size = 0.1, ncol = 3)) + NoLegend() + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
        p2 <- SpatialPlot(sce, features = nCount, pt.size.factor = opts$ptsize)  + theme(legend.position = 'right')
        plot(p1 + p2)
        while (!is.null(dev.list()))  dev.off() 
     
        # 04 plot nFeature_Spatial with SpatialPlot
        pdf(file = paste(opts$out,'/', opts$prefix,'_nFeature_Spatial',  '.pdf', sep = ''), w = 10, h = 5)
        p1  <- suppressWarnings(VlnPlot(sce, group.by = 'orig.ident', features = nFeature, pt.size = 0.1, ncol = 3)) + NoLegend() + theme(axis.text.x=element_blank(), axis.title.x=element_blank())
        p2 <- SpatialPlot(sce, features = nFeature, pt.size.factor = opts$ptsize)  + theme(legend.position = 'right')
        plot(p1 + p2)
        # while (!is.null(dev.list()))  dev.off() 
        dev.off()
    }
}


#' filteration 
filter_sce <- function(sce, opts){
    # 官方教程仅对基因数以及线粒体基因的UMI比例做了限制，其它的参数后续可以自行添加
    # assay ('Spatial', 'RNA')
    if (opts$assays == 'Spatial'){
        selected_c <- WhichCells(sce, expression = (nCount_Spatial >opts$minCount & nCount_Spatial < opts$maxCount & nFeature_Spatial > opts$minFeature & nFeature_Spatial < opts$maxFeature & percent.mt < opts$mtRatio & percent.ribo < opts$rbRatio & percent.hb < opts$hbRatio & percent.plat < opts$platRatio))
    } else if(opts$assays == 'RNA') {
        selected_c <- WhichCells(sce, expression = (nCount_RNA >opts$minCount & nCount_RNA < opts$maxCount & nFeature_RNA > opts$minFeature & nFeature_RNA < opts$maxFeature & percent.mt < opts$mtRatio & percent.ribo < opts$rbRatio & percent.hb < opts$hbRatio & percent.plat < opts$platRatio))
    }
  
    selected_f <- rownames(sce)[Matrix::rowSums(sce) > opts$minCell]
    sce <- subset(sce, features = selected_f, cells = selected_c)
}


#' Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
compute_most_expressed_gene <- function(sce, opts){
    if (opts$assays == 'Spatial'){
      C <- sce@assays$Spatial@counts
    } else if(opts$assays == 'RNA') {
      C <- sce@assays$RNA@counts
    }
    
    C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
    most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
    
    pdf(file = paste(opts$out, '/', opts$prefix, "_most_expressed_gene_boxplot.pdf", sep = ""))
    par(mar = c(4, 8, 2, 1))
    boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
            col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
    while (!is.null(dev.list()))  dev.off() 
}


calc_qc_metrics <- function(sce, opts){
    ## calculate median, mean, std, min 25%, 50%, 75%, max,
    library(psych)
    qc_metrics_data <- data.frame(
                                  percent_mt = c(sce$percent.mt),
                                  percent_rb = c(sce$percent.ribo),
                                  percent_hb = c(sce$percent.hb),
                                  percent_plat = c(sce$percent.plat))
    
    
    if (opts$assays == 'Spatial'){
      qc_metrics_data$umi_counts = c(sce$nCount_Spatial)
      qc_metrics_data$features_num = c(sce$nFeature_Spatial)
    } else if(opts$assays == 'RNA') {
      qc_metrics_data$umi_counts = c(sce$nCount_RNA)
      qc_metrics_data$features_num = c(sce$nFeature_RNA)
    }
    
    qc_metrics <- as.data.frame(matrix(ncol = 11, nrow = 6))
    rownames(qc_metrics) <- c('umi_counts', 'features_num', 'percent_mt',
                              'percent_rb', 'percent_hb', 'percent_plat')
    
    
    for (i in colnames(qc_metrics_data)){
        qc_metrics[i, ] = unlist(describe(qc_metrics_data[, i], na.rm = T, omit = T, skew = F, fast = T, quant = c(.25, .50, .75)))
    }
    colnames(qc_metrics) <- c('vars', 'number','mean', 'sd', 'min', 'max', 'range', 'se', '0.25', '0.5', '0.75')
    # delete the 'vars' and 'range'
    qc_metrics[ , c('vars', 'range')] <- list(NULL)
    # save the results
    write.table(qc_metrics, file = paste(opts$out,'/', opts$prefix, '_qc_metrics_results.txt', sep = ''), sep = '\t',  col.names = T,row.names = T )
    write.table(qc_metrics_data, file = paste(opts$out,'/', opts$prefix, '_qc_metrics_data.txt', sep = ''), sep = '\t', col.names = T, row.names = T )
}


#' calculate percent of mt, ribo, hb and plat each bin
calc_percenterage_feature_set <- function(sce, opts){
    #' 01 calculate percent of mt, ribo, hb and plat each bin
    sce <- PercentageFeatureSet(sce, opts$mtPattern, col.name = "percent.mt")
    sce <- PercentageFeatureSet(sce, opts$rbPattern, col.name = "percent.ribo")
    sce <- PercentageFeatureSet(sce, opts$hbPattern, col.name = "percent.hb")
    sce <- PercentageFeatureSet(sce, opts$platPattern, col.name = "percent.plat")
    #' 02 compute most expressed gene in seurat object
    # cost a lot of memory, skip this step
    # compute_most_expressed_gene(sce, opts)
    #' 03 preview the seurat object
    preview_seuratObject(sce, opts)
    #' 04 calculate qc metrics
    calc_qc_metrics(sce, opts)
    return(sce)
}


##############Future work####################
#' find doublet (to be continue)
#' 
#' Calculate cell-cycle scores
#' 
#' Sample sex (filter sex genes)




# Main work flow
# 第一种情况：只将gem转成seurat文件，同时完成QC统计。
# 第二种情况：导入seurat文件，给定参数进行标准化和过滤。并输出过滤后的统计结果。
# 第三种情况：将gem转成seurat文件，同时直接进行过滤和标准化（不推荐）
if ((! is.null(opts$input)) & (! is.null(opts$out)) ){
    suppressMessages(library(Seurat))
    suppressMessages(library(SeuratObject))
    suppressMessages(library(Matrix))
    # load library for function gem_to_seuratObject()
    suppressMessages(library(data.table))
    suppressMessages(library(rjson))
    suppressMessages(library(ggplot2))

    
    #' 00 create directory 
    dir.create(file.path(opts$out), recursive = TRUE, showWarnings = FALSE)
    setwd(file.path(opts$out))
  
    library(tools)
    # 使用file_ext() 获取文件的扩展名
    suffix <- file_ext(opts$input)
    
    #' 01 create seurat object
    if (!file_test('-f', opts$input)| (opts$fileType == 'matrix' & suffix == 'gz') ) {
        #' i. create seuratObject from a directory. '10X' single cell data format (matrix.mtx, barcodes.tsv, genes.tsv)
        if (!file_test('-f', opts$input)){
        data <- Read10X(opts$input, gene.column = opts$geneColumn)
        #' i. create seuratObject from a *.txt.gz file col: cells, row: genes
        } else if (suffix == 'gz') {
            data <- read.table(gzfile(opts$input), header = T, row.name = 1, check.names = FALSE)
        } 

        sce <- CreateSeuratObject(counts = data, project = opts$prefix,  min.cells = 3, min.features = 200)
        opts$assays <- 'RNA'
        #' ii. calculate percent of mt, ribo, hb and plat each cell
        sce <- calc_percenterage_feature_set(sce, opts)
      
    } else if (suffix == 'gem' | (opts$fileType == 'gem' & suffix == 'gz')){
        #' i. create seuratObject from '*.gem' file
        sce <- gem_to_seuratObject(opts$input, prefix = opts$prefix, binsize = opts$binsize)
        opts$assays <- 'Spatial'
        #' ii. calculate percent of mt, ribo, hb and plat each bin
        sce <- calc_percenterage_feature_set(sce, opts)
      
    } else if (suffix %in% c('rds', 'RDS')){
        sce <- readRDS(opts$input)
        if ('Spatial' == DefaultAssay(sce)|'RNA' == DefaultAssay(sce)){
          opts$mFilter <- 'TRUE'
        } else if ('SCT' == DefaultAssay(sce)){
          stop("The loaded Seurat object already contains SCT assay", call.=FALSE)
        }
        # 如果导入的seurat对象没有做QC的统计，那么先给出线粒体基因等百分比。
        #' 01 calculate percent of mt, ribo, hb and plat each bin
        sce <- PercentageFeatureSet(sce, opts$mtPattern, col.name = "percent.mt")
        sce <- PercentageFeatureSet(sce, opts$rbPattern, col.name = "percent.ribo")
        sce <- PercentageFeatureSet(sce, opts$hbPattern, col.name = "percent.hb")
        sce <- PercentageFeatureSet(sce, opts$platPattern, col.name = "percent.plat")
    } else {
        stop("Please provide a proper file name, e.g. *.gem, *.gem.gz, *.matrix.txt.gz, *.rds, a directory including 10X data", call.=FALSE)
    } 

    
    
    #' 02 filteration and normalization
    if (opts$mFilter) {
        opts$assays <- DefaultAssay(sce)
        sce <- filter_sce(sce, opts)
        ## 仅保留过滤后的原始数据
        saveRDS(sce, paste(opts$out, '/', opts$prefix, '_filt_seuratObject.rds', sep = ''))


        new_opts <- opts
        new_opts$prefix <- paste(opts$prefix, 'filt', sep='_')
        sce <- calc_percenterage_feature_set(sce, new_opts)
        
        # preview_seuratObject(sce, new_opts)
        # calc_qc_metrics(sce, new_opts)
        
        #' 04 normalization with SCTransform
        sce <- SCTransform(sce, assay = opts$assays, new.assay.name = "SCT", verbose = FALSE)

        #' Extract the filtered_norm_matrix and filtered_matrix
        # if (opts$assays == 'Spatial'){
        #   f_out <- gzfile(paste0(opts$out, '/', opts$prefix, '_bin', opts$binsize, '_filt_norm_matrix.txt.gz'), 'w')
        # } else if (opts$assays == 'RNA'){
        #   f_out <- gzfile(paste0(opts$out, '/', opts$prefix, '_filt_norm_matrix.txt.gz'), 'w')
        # }
        if (opts$assays == 'Spatial'){
            f_out_filtered <- gzfile(paste0(opts$out, '/', opts$prefix, '_bin', opts$binsize, '_filt_matrix.txt.gz'), 'w')
            write.table(as.data.frame(sce@assays$Spatial@counts), f_out_filtered, sep = '\t', quote = FALSE)
            close(f_out_filtered)
            f_out <- gzfile(paste0(opts$out, '/', opts$prefix, '_bin', opts$binsize, '_filt_norm_matrix.txt.gz'), 'w')
        } else if (opts$assays == 'RNA'){
            f_out_filtered <- gzfile(paste0(opts$out, '/', opts$prefix, '_filt_matrix.txt.gz'), 'w')
            write.table(as.data.frame(sce@assays$RNA@counts), f_out_filtered, sep = '\t', quote = FALSE)
            close(f_out_filtered)
            f_out <- gzfile(paste0(opts$out, '/', opts$prefix, '_filt_norm_matrix.txt.gz'), 'w')
        }
        write.table(as.data.frame(sce@assays$SCT@counts), f_out, sep = '\t', quote = FALSE)
        close(f_out)

        #' Save the sce rds for further analysis
        # check the variable size
        format(object.size(sce), units = 'auto')
        # set the assay 'Spatial' as NULL to save memory 
        assays_selected <- intersect(names(sce@assays), c('Spatial', 'RNA'))
        # sce[[assays_selected]] <- NULL
        
        ## Display size of seurat object 
        format(object.size(sce), units = 'auto')
        saveRDS(sce, paste(opts$out, '/', opts$prefix, '_filt_norm_seuratObject.rds', sep = ''))

    } else {
        #' Save the sce rds for further analysis
        saveRDS(sce, paste(opts$out, '/', opts$prefix, '_raw_seuratObject.rds', sep = ''))
    }

} else {
    print(parser$print_help())
    stop("You need to provide the following parameters: --input, --out", call.=FALSE)
}
