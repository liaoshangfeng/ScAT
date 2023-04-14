suppressMessages(library(Seurat))
suppressMessages(library(SeuratObject))
suppressMessages(library(SeuratDisk))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
suppressMessages(library(rjson))
suppressMessages(library(RColorBrewer))

###使用时，修改inputfile和prefix即可，自动输出在当前文件夹内
inputfile <- "C5.5w.gem"
prefix='sample'

gem_to_seuratObject <- function(data, prefix = 'sample', binsize = 1){
    #' group counts into bins
    data$x <- trunc(data$x / binsize) * binsize
    data$y <- trunc(data$y / binsize) * binsize
    #trunc取整函数，只取整数部分
    #将横纵坐标根据binsize取整，比如，如果binszie=10,x=2401,y=4522,
    #那么处理后，x=2400,y=4520
    
    if ('MIDCounts' %in% colnames(data)) {
      # aa %in% bb，判断aa是否存在于bb中
      data <- data[, .(counts=sum(MIDCounts)), by = .(geneID, x, y)]
    } else {
      data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
    }
    
    #' create sparse matrix from stereo
    data$cell <- paste0(prefix, ':', data$x, '_', data$y)
    data$geneIdx <- match(data$geneID, unique(data$geneID))
    data$cellIdx <- match(data$cell, unique(data$cell))
    #match(aa,bb)返回bb在aa中的位置
    
    # if (! is.null(opts$binsize)){
    #   #is.null(aa) 判断aa是否为空，如果空，返回Ture；
    #   #! is.null(aa)，相反地，判断aa是否非空，如果非空，返回True；
    #   write.table(data, file = paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '.tsv'), 
    #               quote = FALSE, sep = '\t', row.names = FALSE)
    # }
    #write.table(aa,file=xx,sep =" ",quote=TRUE,row.names =TRUE, col.names =TRUE)
    #把内容aa输出到文件xx中
    
    mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                        dimnames = list(unique(data$geneID), unique(data$cell)))
    #sparseMatrix稀疏矩阵函数
    
    cell_coords <- unique(data[, c('cell', 'x', 'y')])
    #unique去除重复函数，删除cell,x和y都一样的行
    
    rownames(cell_coords) <- cell_coords$cell
    
    #cell_coords$cell <- NULL
    
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
    
    spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                      scale.factors = fromJSON(scalefactors_json), 
                                      tissue.positions = tissue_positions_list)
    #可以理解为构建一个spatial背景
    
    #' import image into seurat object
    spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
    DefaultAssay(spatialObj) <- 'Spatial'
    
    seurat_spatialObj[['slice1']] <- spatialObj
    rm("spatialObj")
    rm("data")
    rm("mat")
    
    return(seurat_spatialObj)
    }

#' seurat clustering workflow                                                                                                                                                                                                                  
Clustering <- function(object, dims = 30, resolution = 0.8){                                                                                                                                                                                                     
  object <- RunPCA(object)                                                                                                                                                                                                                   
  object <- FindNeighbors(object, dims = 1:dims)                                                                                                                                                                                             
  object <- FindClusters(object, verbose = FALSE, resolution = resolution)                                                                                                                                                              
  object <- RunUMAP(object, dims = 1:dims)                                                                                                                                                                                                   
  return(object)                                                                                                                                                                                                                             
}

# read the data
#out_dir <- "/home/tonyleao/wkd/data_to_viz/01.preprocessing/"
# file1 <- "DP8400016190TL_C5.bin200.Lasso.gem"
# file2 <- "FP200000477TR_E4.bin200.Lasso.gem"
# file3 <- "FP200000477TR_F4.bin200.Lasso.gem"
# file4 <- "DP8400013846TR_F5_bin200.gem"
#file1 <- "C5.5w.gem"


data <- fread(file = file1)

ob1 <- gem_to_seuratObject(data, binsize = 1, prefix = prefix_name)
saveRDS(ob1, paste(prefix_name,'.seuratObject.rds')
