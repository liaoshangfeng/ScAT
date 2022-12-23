
library('Seurat')
library('SPARK')
library('Matrix')
library('optparse')
library('ggplot2')


option_list = list(
                   make_option(c("-i", "--input"), type="character", default=NULL, dest = "input",
                               help="Gene expression count matrix", metavar="character"),
                   make_option(c("-o", "--out"), type="character", default=NULL, dest = "out",
                               help="result file", metavar="character"),
                   make_option(c("-g", "--genefile"), type="character", default=NULL, dest = "genefile",
                   help="gene file", metavar="character")
                   )

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
#rawcount=read.table(opt$input,check.names = FALSE)
Stdata=readRDS(opt$input)
rawcount=as.data.frame(Stdata@assays$Spatial@counts)
rawcount = as.matrix(rawcount)
rawcount <- as(rawcount, "dgTMatrix")
genefile=opt$genefile
outdir=opt$out


  if (!is.null(outdir)){
    if (! file.exists(opt$out)){
        dir.create(opt$out, recursive = TRUE)}
        }


##获取坐标信息
Get_coord_info<-function(){
  coord_col=sapply(strsplit(colnames(rawcount),split=":"),"[",2)
  info <- cbind.data.frame(x=as.numeric(sapply(strsplit(coord_col,split="_"),"[",1)),
                           y=as.numeric(sapply(strsplit(coord_col,split="_"),"[",2)))
  rownames(info)  <- colnames(rawcount)
  location  <- as.matrix(info) 
  return(location)
}


spark_calculate<-function(){
      

      location<-Get_coord_info()
      
      mt_idx  <- grep("mt-",rownames(rawcount))
                 if(length(mt_idx)!=0){
      rawcount<- rawcount[-mt_idx,]
        }
     
       sparkX <- sparkx(rawcount,location,numCores=1,option="mixture")
      }

if  (!is.null(genefile)){
     
     cluster_name_list <- function(genefile){
      data=read.table(genefile)
       for (i in seq(dim(data)[1])){
            if (i==1){
            cluster_name= data[i,]
       }else{
        cluster_name=paste0(cluster_name,',',data[i,])}
        }
        cluster_name <- unlist(strsplit(cluster_name,split=","))
        print(cluster_name)
        return (cluster_name)
        }
 
     p=Seurat::SpatialPlot(Stdata,features=cluster_name_list(genefile),
        label = TRUE, pt.size.factor=4,alpha = c(0.3, 1)) 
     ggsave(p,file=file.path(opt$out,'spatial.pattern.gene.pdf'),width=10,height = 8)
     } 
     
  if (is.null(genefile)) {
    spark_res=spark_calculate()
    spark_res=as.data.frame(spark_res)
    spark_res=spark_res[order(spark_res[,'res_mtest.adjustedPval']),]
    write.table(spark_res,file=file.path(opt$out,'spatial.pattern.gene.tsv'),quote=F,sep='\t') 
    }  
       



