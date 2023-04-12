suppressPackageStartupMessages({
  library(optparse)
  library(SingleCellExperiment)
  library(ggplot2)
  library(BayesSpace)
  library(Seurat)
})
set.seed(100)


## Parameters
parser = argparse::ArgumentParser(description = 'Spatial transcriptome cluster with BayesSpace')
parser$add_argument('--Stdata', dest = 'Stdata', default=NULL, help = 'Input spatial transcriptome data')
parser$add_argument('--clusterNumber', dest = 'clusterNumber',  default=NULL, help = 'Set cluster number of BayesSpace spatial cluster')
parser$add_argument('--iterations', dest = 'iterations',  default=NULL, help = 'iterations number of MCMC')
parser$add_argument('--out', dest = 'out',default=NULL, help = 'Output file into this directory')
opts = parser$parse_args()


## Functions

#' define a function for data preparation
input_format <- function(){
    #rds_file='/jdfssz1/ST_TSCBI/PROJECT_temp/USER/huangke/10.STspotlight.refine/deconv.v4/ST.ms.gem.raw.coord.rds'
    rds_file=opts$Stdata
    Stdata=readRDS(rds_file)
    df=data.frame(X=Stdata@meta.data$x,Y=Stdata@meta.data$y)
    df$spot_id=paste(df$X, df$Y, sep="X")
    df<-df[, c("spot_id","X","Y")]
    write.table(df,file=file.path(opts$out,'00.bin.info.tsv'),quote = FALSE,col.names = FALSE,row.names = FALSE,sep='\t')
    coord_file=file.path(opts$out,'00.bin.info.tsv')
    counts = as(Stdata@assays$Spatial@counts, "dgCMatrix")

    colData = read.table( coord_file, stringsAsFactors=F, sep="\t" )
    names(colData) = c( "spot", "row", "col" )
    rownames(colData) = colData$spot
    genes = rownames(counts)
    rowData = data.frame(genes)

    colnames(counts) <- colData$spot
    return( list( one=counts, two=rowData,three=colData ) )
}

#' define a function for spatial cluster with Bayes
bayes_spatial_cluster <- function(){
      fun_res=input_format()
      counts=fun_res$one
      rowData=fun_res$two
      colData=fun_res$three
      ##构建SingleCellExperiment objects
      sce <- SingleCellExperiment(assays=list(counts=counts),
                                rowData=rowData,
                                colData=colData)

       
       ##数据预准备，包括log标准化和PCA
      sce <- spatialPreprocess(sce, platform="ST", 
                               n.PCs=30, n.HVGs=3000, log.normalize=TRUE)
       
      ##对spot进行聚类，并添加聚类标签到 SingleCellExperiment 对象
      ##q：指定聚类的数目
      ##d：使用的top PCA数
      ##init.method 初始聚类方法 mclust or kmeans
      ## model t or normal
      ##gamm 平滑参数 (1-3即可)
      ##nrep 隐马迭代次数
      ## burn.in 相当于去除前200次初始迭代
       
      sce <- spatialCluster(sce, q=as.numeric(opts$clusterNumber), platform="ST", d=30,
                           init.method="mclust", model="t", gamma=2,
                           nrep=as.numeric(opts$iterations), burn.in=200,
                           save.chain=FALSE)
  
    out.df = colData(sce)
    write.table(out.df, file=file.path(opts$out,'02.cluster.info.txt'), quote=F, row.names=T, col.names=T, sep="\t")

}


#' define a function for color palette
ColorPalette <- function(number) {
  print(number)
  if (number < 25) {
    colorScheme = c("#999999","#FF0099","#E69F00","#56B4E9","#009E73","#F0E442",
                    "#0072B2","#D55E00","#CC79A7","#990000","#9900cc","#66FF66",
                    "#663300","#0000FF","#CC0033","#FF0000","#000099","#660066",
                    "#333333","#FFCCCC","#993366","#33CC33","#000099","#CC9900")
  } else {
    suppressMessages(library(RColorBrewer))
    qual_col_pals <-
      brewer.pal.info[brewer.pal.info$category == 'qual', ]
    colorScheme <-
      unique(unlist(mapply(
        brewer.pal,
        qual_col_pals$maxcolors,
        rownames(qual_col_pals)
      )))
    return (colorScheme)
  }
}

#' define a function to visualize the spatial cluster
bayes_cluster_figure <- function(point_size = 1){
  df <- read.table(file.path(opts$out, '02.cluster.info.txt'), header = T, row.names = 1)
  colnames(df) <- c('bin_id', 'x', 'y', 'B', 'cluster.init', 'spatial.cluster')
  df$y <- -(df$y)  ## you can comment this for your need
  df$spatial.cluster <- as.character(df$spatial.cluster)
  colorPalette <- ColorPalette(length(unique(df$spatial.cluster)))
  
  p0 <-
    ggplot(df, aes(x = x, y = y, color = spatial.cluster)) +
    geom_point(shape = 19, size = point_size) +  ## you can modify the point size to avoid overlap
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = 'right'
    ) +
    coord_fixed()  +
    scale_color_manual(values = colorPalette) +
    guides(colour = guide_legend(
      override.aes = list(size = 3),
      nrow = 15,
      title = 'Clusters'
    ))  +
    theme_void()   
  
  pdf(file = file.path(opts$out,'03.BayesSpace.spatial.cluster.pdf'), width = 6, height = 6)
  p0
  dev.off()
  
}


if (!is.null(opts$out)){   
    if (! file.exists(opts$out)){ 
        dir.create(opts$out, recursive = TRUE)} 
}
##调用函数 
bayes_spatial_cluster()
bayes_cluster_figure()

