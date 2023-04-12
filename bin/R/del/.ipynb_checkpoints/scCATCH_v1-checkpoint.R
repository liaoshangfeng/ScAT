##脚本功能:用scCATCH软件对分类结果进行注释

suppressMessages(library(optparse))
suppressMessages(library(scCATCH))
suppressMessages(library(Seurat))

option_list = list(
                make_option(c("-i", "--inputrds"), type="character", default=NULL, dest = "input",
                               help="input dirname", metavar="character"), 
                make_option(c("-s", "--species"), type="character", default="mouse", dest = "species",
                               help="species name (Human or Mouse)", metavar="character"),
                make_option(c("-t", "--tissue"), type="character", default=NULL, dest = "tissue",
                               help="tissue name", metavar="character"),
                make_option(c("-c", "--cancer"), type="character", default=NULL, dest = "cancer",
                               help="cancer name", metavar="character"),
                make_option(c("--cluster"), type="character", default="All", dest = "cluster",
                               help="choose the clusters that you want to annotate,eg:c('1','2')", metavar="character"),
                make_option(c("--cell_min_pct"), type="numeric", default = 0.25),
                make_option(c("--logfc"), type="numeric", default = 0.25),
                make_option(c("--pvalue"), type="numeric", default = 0.05),
                make_option(c("-o", "--out"), type="character", default=NULL, dest = "out",
                               help="directory to save file [default= %default]", metavar = "character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


##函数功能：对聚类后的数据进行注释
sc_anno <- function(rds, species, tissue, cancer, cluster, cell_min_pct, logfc, pvalue){
    object <- readRDS(rds)
        # change the Assays to 'RNA'
        if (DefaultAssay(object)=='Spatial'){
            DefaultAssay(object)
            object <- RenameAssays(object, Spatial='RNA')
        }

    #根据表达量找到marker genes
    clu_markers <- suppressMessages(findmarkergene(object = object, 
	species = species, 
	tissue = tissue, 
	cancer = cancer, 
	cluster = cluster, 
	match_CellMatch = TRUE, 
	cell_min_pct = cell_min_pct, 
	logfc = logfc, 
	pvalue = pvalue))
    #结合高变基因和文献信息对细胞簇进行注释
    clu_ann <- suppressMessages(scCATCH(object = clu_markers$clu_markers, species = species, cancer = cancer, tissue = tissue))
    #对细胞簇进行重命名注释
    new.cluster.ids <- clu_ann$cell_type
    names(new.cluster.ids) <- clu_ann$cluster
    object <- RenameIdents(object, new.cluster.ids)
    return(list(anno=clu_ann, object=object))
}


MainPlotData<-function(){

    #00.对输入的rds文件进行注释
    result <- sc_anno(rds = opt$input, 
	species = opt$species, 
	cancer = opt$cancer,
	tissue = opt$tissue, 
	cluster = opt$cluster, 
	cell_min_pct = opt$cell_min_pct, 
	logfc = opt$logfc, 
	pvalue = opt$pvalue)
    anno <- result$anno
    sce <- result$object
    #01.保存注释信息：A data.frame containing matched cell type for each cluster, related marker genes, evidence based score and PMID
    write.table(anno, file = paste(opt$out, "cluster_anno.txt", sep = "/"), sep = "\t")
    
    #02.输出注释后的聚类图
    pdf(file = paste(opt$out, "anno_cluster.pdf", sep = "/"))
    plot <- DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    print (plot)
    while (!is.null(dev.list()))  dev.off()
    
    #03.保存注释后的rds的文件
    saveRDS(sce, paste(opt$out, "seuratObject_anno.rds", sep = "/")) 
}

##新建输出文件夹
if (! file.exists(opt$out)){
    dir.create(opt$out)
}
##调用主函数
scCATCH_anno=MainPlotData()
