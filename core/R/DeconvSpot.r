library(optparse)
library(Matrix)
library(data.table)
library(Seurat)
library(dplyr)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(ggplot2)
suppressMessages(require(imager))
suppressMessages(require(png))
suppressMessages(require(jpeg))
suppressMessages(require(grid))




parser = argparse::ArgumentParser(description = 'Deconvolute spatial transcriptomics spots with single-cell transcriptomes')
parser$add_argument('--Scdata', dest = 'Scdata', default=NULL, help = 'Input normalized data (*.txt.gz)')
parser$add_argument('--Stdata', dest = 'Stdata', default=NULL, help = 'Input spatial transcriptome data')
#parser$add_argument('--StHE', dest = 'StHE', default=NULL,help = 'HE png of spatial transcriptome')
parser$add_argument('--ClusterFile', dest = 'clusterFile',  default=NULL, help = 'View a set of specific cluster locations')
parser$add_argument('--piescale', dest = 'piescale',default=0.8, help = 'Pie size of deconv plot')
parser$add_argument('--outdir', dest = 'outdir',default=NULL, help = 'Output file into this directory')

opts = parser$parse_args()
Scdata<- readRDS(opts$Scdata)
Stdata<-readRDS(opts$Stdata)
clusterFile=opts$clusterFile
piescale= as.numeric(opts$piescale)
img_path='/jdfssz1/ST_TSCBI/PROJECT_temp/USER/huangke/13.brain.project/09.deconvSpot/01.deconv.E4/blank.size.png'
outdir=opts$outdir 


#' 01 格式转换;clusterFile内容转成向量，用于传入函数
cluster_name_list <- function(clusterFile){
      data=read.table(clusterFile)
       for (i in 1:dim(data)[1]){
            if (i==1){
            cluster_name= data$V1[i]
       }
      else{
        cluster_name=paste0(cluster_name,',',data$V1[i])}
        }
        cluster_name <- unlist(strsplit(cluster_name,split=","))
     
        return (cluster_name)
        }


#' 03  函数用来根据单细胞和空间数据做umap图

Scdata_rds <- function(Scdata) {
  Scdata <- Seurat::SCTransform(Scdata, verbose = FALSE) %>%
  Seurat::RunPCA(verbose = FALSE) %>%
  Seurat::FindNeighbors(dims = 1:30, verbose = FALSE)%>%
  Seurat::FindClusters(verbose = FALSE)%>%
  Seurat::RunUMAP(dims = 1:30, verbose = FALSE)
  p0<-Seurat::DimPlot(Scdata,
          label = TRUE) + Seurat::NoLegend()
  ggsave(p0,file=file.path(outdir,'01.single.cell.cluster.pdf'),width=10,height=8)
  return(Scdata)
  }

Stdata_rds <- function(Stdata) {
  Stdata <- Seurat::SCTransform(Stdata, assay = "Spatial", verbose = FALSE) %>% 
  Seurat::RunPCA(assay='SCT',verbose = FALSE) %>%
  Seurat::FindNeighbors(dims = 1:30, verbose = FALSE)%>%
  Seurat::FindClusters(verbose = FALSE)%>%
  Seurat::RunUMAP(dims = 1:30, verbose = FALSE) 
  p0<-Seurat::DimPlot(Stdata,
                      label = TRUE) + labs(fill=guide_legend(title="New Legend Title"))
  p1<-Seurat::SpatialPlot(Stdata,
                      label = TRUE, pt.size.factor=4,alpha = c(0.3, 1)) + Seurat::NoLegend()
  ggsave(p0+p1,file=file.path(outdir,'02.Spatial.cluster.pdf'),width=12,height=8)
  
  return(Stdata)
}


#‘ 06 解卷积主函数:根据单细胞cell marker在每个细胞类型中的topic，解析spot细胞组成类型
Deconv_Spot_fun<-function(Scdata,Stdata){
    #’ 07  寻找单细胞数据的cell marker
    #Seurat::Idents(object = Scdata) <- Scdata@meta.data$subclass
    Seurat::Idents(object = Scdata) <- Scdata@meta.data$seurat_clusters
    cluster_markers_all <- Seurat::FindAllMarkers(object = Scdata, 
                                              assay = "RNA",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)

    saveRDS(object = cluster_markers_all,file.path(outdir,'03.ScDATA.markers.rds'))
    #cluster_markers_all=readRDS(file.path(outdir,'03.ScDATA.markers.rds'))
    
    #‘ 08 根据单细胞数据，单细胞cell marker数据，推断空间数据每个spot的细胞组成
    set.seed(123)
    spotlight_ls <- spotlight_deconvolution(
      se_sc = Scdata,
      counts_spatial = Stdata@assays$Spatial@counts,
      clust_vr = "seurat_clusters", # Variable in sc_seu containing the cell-type annotation
      cluster_markers = cluster_markers_all, # Dataframe with the marker genes
      cl_n = 100, # number of cells per cell type to use
      hvg = 2500, # Number of HVG to use
      ntop = NULL, # How many of the marker genes to use (by default all)
      transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
      method = "nsNMF", # Factorization method
      min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
      )
      #' 存储解卷积结果
      saveRDS(object = spotlight_ls, file = file.path(outdir,"04.spotlight_ls.rds"))
      }

##09 添加解卷积结果信息到stdata

St_spot_plot_format<-function(){
  spotlight_ls=readRDS(file.path(outdir,"04.spotlight_ls.rds"))
  nmf_mod <- spotlight_ls[[1]]
  decon_mtrx <- spotlight_ls[[2]]
  decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
  decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
  rownames(decon_mtrx) <- colnames(Stdata)
  
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  
  Stdata@meta.data <- Stdata@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  return(Stdata)
}

## 10 获取细胞类型列表
#Cell_type_all <- function(Stdata){
#  spotlight_ls=readRDS(file.path(outdir,"spotlight_ls.rds"))
#  decon_mtrx <- spotlight_ls[[2]]
#  decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
#  decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
#  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"]) 
  
#  cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
#  return(cell_types_all)
#}

##  Plot  figure  acoording to deconv result   
Deconv_Scatterpie<-function(se_obj,cell_types_all,cell_types_interest,pie_scale){
     scatterpie_alpha = 1
     pie_scale = pie_scale
    
    # Check variables
    if (!is(se_obj, "Seurat")) stop("ERROR: se_obj must be a Seurat object!")
    if (! is(cell_types_all, "vector")) stop("ERROR: cell_types_all must be a vector/list object!")
    if (!is.character(img_path)) stop("ERROR: must be a character string!")
    if (!(is(cell_types_interest, "vector") | is.null(cell_types_interest))) stop("ERROR: cell_types_interest must be a vector/list object or NULL!")
    if (!is.numeric(scatterpie_alpha)) stop("ERROR: scatterpie_alpha must be numeric between 0 and 1!")
    if (!is.numeric(pie_scale)) stop("ERROR: pie_scale must be numeric between 0 and 1!")
    
    metadata_ds <- data.frame(se_obj@meta.data)
    colnames(metadata_ds) <- colnames(se_obj@meta.data)
    
    if (is.null(cell_types_interest)) {
      cell_types_interest <- cell_types_all
    }
    
    # If not all cell types are in the cell types of interest we only want to keep those spots which have at least one of the cell types of interest
    if (!all(cell_types_all %in% cell_types_interest)) {
      print(cell_types_interest)
      metadata_ds <- metadata_ds %>%
        tibble::rownames_to_column("barcodeID") %>%
        dplyr::mutate(rsum = base::rowSums(.[, cell_types_interest,
                                             drop = FALSE])) %>%
        dplyr::filter(rsum != 0) %>%
        dplyr::select("barcodeID") %>%
        dplyr::left_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                         by = "barcodeID") %>%
        tibble::column_to_rownames("barcodeID")
    }
    
    ## If slice is not selected set it to the first element in the list of slices
  #  if (is.null(slice) | (!is.null(slice) && !slice %in% names(se_obj@images))) {
      slice <- names(se_obj@images)[1]
   #   warning(sprintf("Using slice %s", slice))
   # }
    
    ## Preprocess data
    spatial_coord <- data.frame(se_obj@images$slice@coordinates) %>%
      tibble::rownames_to_column("barcodeID") %>%
      dplyr::mutate(
        imagerow_scaled =
          imagerow * se_obj@images$slice@scale.factors$lowres,
        imagecol_scaled =
          imagecol * se_obj@images$slice@scale.factors$lowres) %>%
      dplyr::inner_join(metadata_ds %>% tibble::rownames_to_column("barcodeID"),
                        by = "barcodeID")
    
    ### Load histological image into R
    #### Extract file format, JPEG or PNG
    img_frmt <- base::tolower(stringr::str_sub(img_path, -4, -1))
    
    if(img_frmt %in% c(".jpg", "jpeg")) {
      img <- jpeg::readJPEG(img_path)
    } else if (img_frmt == ".png") {
      img <- png::readPNG(img_path)
    }
    
    # Convert image to grob object
    img_grob <- grid::rasterGrob(img,
                                 interpolate = FALSE,
                                 width = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc"))
    
    
    #if (cell_types_all[1]=='0'){
    #  cell_type_col=paste("X", (1:(length(cell_types_all))-1), sep = "")
    #}else{
    #  cell_type_col=cell_types_all
    #}
    # color_list <- c('#FB9A99','darkturquoise','#FDBF6F','deeppink1','orchid1','#A52A2A', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
    #                 'skyblue2', 'palegreen2', '#CAB2D6', 'gray70', 'khaki2', 
    #                 'maroon', 'steelblue4', 'green1', 'yellow4', 'yellow3','darkorange4','dodgerblue2')

    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
    color_list <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))) 
    all_col=length(colnames(spatial_coord))
    colnames(spatial_coord)[(all_col-length(cell_types_all)):(all_col-1)]=cell_types_all

    #colnames(spatial_coord)[19:(19+length(cell_types_all)-1)]=cell_types_all
    write.table(spatial_coord,file=file.path(outdir,"spatial_coord.prop.pip.res.tsv"),quote=F,row.names = F,sep='\t')
    print(spatial_coord)
    ## Plot spatial scatterpie plot
    scatterpie_plt <- suppressMessages(
      ggplot2::ggplot() +
        ggplot2::annotation_custom(
          grob = img_grob,
          xmin = 0,
          xmax = ncol(img),
          ymin = 0,
          ymax = -nrow(img)) +
        scatterpie::geom_scatterpie(
          data = spatial_coord,
          ggplot2::aes(x = imagecol_scaled,
                       y = imagerow_scaled),
          cols = cell_types_all ,
          color = NA,
          alpha = scatterpie_alpha,
          pie_scale = pie_scale)+scale_fill_manual(values=color_list[1:length(cell_types_all)])+
        ggplot2::scale_y_reverse() +
        #ggplot2::ylim(nrow(img), 0) +
        #ggplot2::xlim(0, ncol(img)) +
        cowplot::theme_half_open(11, rel_small = 1) +
        ggplot2::theme_void() +
        ggplot2::coord_fixed(ratio = 1,
                             xlim = NULL,
                             ylim = NULL,
                             expand = TRUE,
                             clip = "on"))
    return(scatterpie_plt)
  }



Main_Fun <- function(){

  #'  若目录下没有输出目录文件夹,新建输出文件夹
  if (!is.null(outdir)){
    if (! file.exists(opts$out)){
        dir.create(opts$out, recursive = TRUE)}
        }

if  (is.null(clusterFile)){

  Scdata<- readRDS(opts$Scdata)
  Stdata<-readRDS(opts$Stdata)
  if (is.null(Stdata@meta.data$seurat_clusters)){
    Stdata=Stdata_rds(Stdata)
    }
       
  #clusterFile=opts$clusterFile
  #piescale=opts$piescale
 # img_path='/jdfssz1/ST_TSCBI/PROJECT_temp/USER/huangke/13.brain.project/09.deconvSpot/01.deconv.E4/blank.size.png'
 
 ##调用解卷积函数进行解卷积
  if (is.null(Scdata@meta.data$seurat_clusters)){
      
         Scdata=Scdata_rds(Scdata)
  }

  Deconv_Spot_fun(Scdata,Stdata)
  
  ##根据解卷积结果文件，将画图信息整理进stdata中，以及获取簇名
  spotlight_ls=readRDS(file.path(outdir,'04.spotlight_ls.rds'))
  nmf_mod <- spotlight_ls[[1]]
  decon_mtrx <- spotlight_ls[[2]]
  decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
  decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
  rownames(decon_mtrx) <- colnames(Stdata)
  
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  
  Stdata@meta.data <- Stdata@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"]) 
  cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
  print(cell_types_all)
  ## 所有bin的解卷积后，细胞类型成分饼图
  p1=Deconv_Scatterpie(Stdata,cell_types_all,NULL,piescale)  
  ggsave(p1,file=file.path(outdir,'05.StData.All.cluster.pie.pdf'),width=10,height=8)
  #Stdata=Stdata_rds(Stdata)
  }
  ##如果指定特定cluster
  if  (!is.null(clusterFile)){
         print('cluster fun')
         clusterNames=cluster_name_list(clusterFile)
         spotlight_ls=readRDS(file.path(outdir,'04.spotlight_ls.rds')) 
         if (is.null(Stdata@meta.data$seurat_clusters)){ 
                 Stdata=Stdata_rds(Stdata)
                 }
        nmf_mod <- spotlight_ls[[1]]
       decon_mtrx <- spotlight_ls[[2]]
      decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
      decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
       decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
       rownames(decon_mtrx) <- colnames(Stdata)
  
       decon_df <- decon_mtrx %>%
            data.frame() %>%
       tibble::rownames_to_column("barcodes")
  
      Stdata@meta.data <- Stdata@meta.data %>%
         tibble::rownames_to_column("barcodes") %>%
         dplyr::left_join(decon_df, by = "barcodes") %>%
         tibble::column_to_rownames("barcodes")
       decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"]) 
        cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]

       #' 12 画出指定细胞类型在空间spot中的比例热图
        p3<-Seurat::SpatialFeaturePlot(
        object = Stdata,
        features = clusterNames,
        image.alpha = 0.1,alpha = c(0, 0.8),
        pt.size.factor=6,ncol=4)&scale_fill_continuous(low="dodgerblue2",high="#A52A2A",limits = c(0,1))
        ggsave(p3,file=file.path(outdir,'06.StData.specific.cluster.prop.pdf'),width=15,height=12)
  
       #’ 13 画出指定细胞类型在空间spot中的比例饼图
       #p4=Deconv_Scatterpie(Stdata,cell_types_all,clusterNames,0.08)
       #ggsave(p4,file=file.path(outdir,'07.StData.specific.cluster.pie.pdf'),width=10,height=8)
}
}


Main_Fun()
