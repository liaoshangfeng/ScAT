f## Cell type annotation with scCATCH

### Get the parameters
parser = argparse::ArgumentParser(description = 'Perform cell type annotation by Seurat')
parser$add_argument('-i', '--inputrds', dest = 'inputrds', help = 'Input raw_seuratObj.rds')
parser$add_argument('-s', '--species', dest = 'species', default='Human', help = 'Species of data (limit "Human", "Mouse") [default %(default)s]')
parser$add_argument('-t', '--tissue', dest = 'tissue', default='NULL', help = 'tissue name [default %(default)s]')
parser$add_argument('-c', '--cancer', dest = 'cancer', default='NULL', help = 'cancer name [default %(default)s]')

parser$add_argument('--cluster', dest = 'cluster', default='All', help = 'choose the clusters that you want to annotate,eg:c("1","2") [default %(default)s]')

parser$add_argument('--cell_min_pct', dest = 'cell_min_pct', default='0.25', type = 'double', help = 'min cell percent [default %(default)s]')
parser$add_argument('--logfc', dest = 'logfc', default='0.25', type = 'double', help = 'logFC value [default %(default)s]')
parser$add_argument('--pvalue', dest = 'pvalue', default='0.25',type = 'double',  help = 'pvalue [default %(default)s]')
parser$add_argument('-o', '--out', dest = 'out', default=NULL, help = 'Directory to save file [default %(default)s]')

opts = parser$parse_args()
print(opts)


#函数功能：对聚类后的数据进行注释
sc_anno <- function(rds, species, tissue,  cluster = 'All', cancer = 'Normal',cell_min_pct, logfc, pvalue){
    object <- readRDS(rds)
    # change the Assays to 'RNA'
    if (DefaultAssay(object)=='Spatial'){
        DefaultAssay(object)
        object <- RenameAssays(object, Spatial='RNA')
    }

    
  object <- scCATCH::createscCATCH(object[["RNA"]]@data, cluster = seurat_cluster) 
  
  #根据表达量找到marker genes
  clu_markers <- suppressMessages(findmarkergene(object = object,
                                                 species = species,
                                                 tissue = tissue,
                                                 cancer = cancer,
                                                 cluster = cluster,
                                                 cell_min_pct = cell_min_pct,
                                                 logfc = logfc,
                                                 pvalue = pvalue))
    
  #结合高变基因和文献信息对细胞簇进行注释
  # clu_ann <- suppressMessages(scCATCH(object = clu_markers$clu_markers, species = species, cancer = cancer, tissue = tissue))
   clu_ann <- scCATCH::findcelltype(clu_markers) 
   scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype, by = "cluster") %>% tibble::column_to_rownames("cell")
   object = AddMetaData(object, scCATCH_ann_data)
    
  #对细胞簇进行重命名注释
  new.cluster.ids <- clu_ann$cell_type
  names(new.cluster.ids) <- clu_ann$cluster
  object <- RenameIdents(object, new.cluster.ids)
  pdf(file = paste(opts$out, "anno_cell_type.pdf", sep = "/"))  
  DimPlot(object, pt.size = 1, label = T, label.size = 4, repel = T, group.by = "cell_type")
  dev.off()
  return(list(anno=clu_ann, object=object))
}

MainPlotData<- function(opts){
    #00.对输入的rds文件进行注释
    result <- sc_anno(rds = opts$input, species = opts$species, cancer = opts$cancer,tissue = opts$tissue, 
                      cluster = opts$cluster, cell_min_pct = opts$cell_min_pct, logfc = opts$logfc, pvalue = opts$pvalue)
    anno <- result$anno
    sce <- result$object
    #01.保存注释信息：A data.frame containing matched cell type for each cluster, related marker genes, evidence based score and PMID
    write.table(anno, file = paste(opts$out, "cluster_anno.txt", sep = "/"), sep = "\t")
    
    #02.输出注释后的聚类图
    pdf(file = paste(opts$out, "anno_cluster.pdf", sep = "/"))
    plot <- DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
    print (plot)
    while (!is.null(dev.list()))  dev.off()
    
    #03.保存注释后的rds的文件
    saveRDS(sce, paste(opts$out, "seuratObject_anno.rds", sep = "/")) 
}

suppressPackageStartupMessages({
    library(scCATCH) 
    library(Seurat)
})
##新建输出文件夹
if (! file.exists(opts$out)){
    dir.create(opts$out)
}
##调用主函数
MainPlotData(opts)