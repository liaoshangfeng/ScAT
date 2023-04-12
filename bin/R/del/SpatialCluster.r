library(optparse)
library(SingleCellExperiment)
library(ggplot2)
library(BayesSpace)
library(Seurat)
set.seed(100)

parser = argparse::ArgumentParser(description = 'Spatial transcriptome cluster with BayesSpace')
parser$add_argument('--Stdata', dest = 'Stdata', default=NULL, help = 'Input spatial transcriptome data')
parser$add_argument('--clusterNumber', dest = 'clusterNumber',  default=NULL, help = 'Set cluster number of BayesSpace spatial cluster')
parser$add_argument('--iterations', dest = 'iterations',  default=NULL, help = 'iterations number of MCMC')
parser$add_argument('--outdir', dest = 'outdir',default=NULL, help = 'Output file into this directory')
opts = parser$parse_args()

input_format<-function(){
    #rds_file='/jdfssz1/ST_TSCBI/PROJECT_temp/USER/huangke/10.STspotlight.refine/deconv.v4/ST.ms.gem.raw.coord.rds'
    rds_file=opts$Stdata
    Stdata=readRDS(rds_file)
    df=data.frame(X=Stdata@meta.data$x,Y=Stdata@meta.data$y)
    df$spot_id=paste(df$X, df$Y, sep="X")
    df<-df[, c("spot_id","X","Y")]
    write.table(df,file=file.path(opts$outdir,'00.bin.info.tsv'),quote = FALSE,col.names = FALSE,row.names = FALSE,sep='\t')
    coord_file=file.path(opts$outdir,'00.bin.info.tsv')
    counts = as(Stdata@assays$Spatial@counts, "dgCMatrix")
    colData = read.table( coord_file, stringsAsFactors=F, sep="\t" )
    names(colData) = c( "spot", "row", "col" )
    rownames(colData) = colData$spot
    genes = rownames(counts)
    rowData = data.frame(genes)
    return( list( one=counts, two=rowData,three=colData ) )
}

bayes_spatial_cluster<-function(){
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
                           save.chain=TRUE)
  


    out.df = colData(sce)
    write.table(out.df, file=file.path(opts$outdir,'02.cluster.info.txt'), quote=F, row.names=T, col.names=F, sep="\t")

}





bayes_cluster_figure<-function(){
    cluster_info <- read.table(file.path(opts$outdir,'02.cluster.info.txt'),header = F)
    cluster_info$X = cluster_info$V3 
    cluster_info$Y = cluster_info$V4
    cluster_info$cluster1 = as.character(cluster_info$V6)

    colorPalette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                  'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                  'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                  'green1', 'yellow4', 'yellow3','darkorange4', 'brow')


     p1<-ggplot(data = cluster_info)+
       geom_point(aes(x=X,y=Y, color = cluster1)) +
      scale_colour_manual(name="",  
                      values = colorPalette)+
      theme_bw() + 
      theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  # theme(legend.position = 'none')
     theme(
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))
    ggsave(p1,file=file.path(opts$outdir,'03.BayesSpace.spatial.cluster.pdf'),width=10,height=8)
}


if (!is.null(opts$outdir)){   
    if (! file.exists(opts$out)){ 
        dir.create(opts$out, recursive = TRUE)} 
}
##调用函数 
bayes_spatial_cluster()
bayes_cluster_figure()

