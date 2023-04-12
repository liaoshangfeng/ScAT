## cell clustering analysis by Seurat for gem data

### Get the parameters
parser = argparse::ArgumentParser(description = 'Script for cell/spatial BIN clustering by Seurat')
parser$add_argument('-i', '--input', dest = 'input', help = 'Input raw_seuratObj.rds')
parser$add_argument('-o', '--out', dest = 'out', help = 'directory where to save the output files, all output files will be indexed by sample ID')

parser$add_argument('--ptsize', dest = 'ptsize', default = 1, type = 'double', help = 'point size for spatial dim plot [default %(default)s]')
parser$add_argument('--prefix', dest = 'prefix', default = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident [default %(default)s]')


parser$add_argument('--pc', dest = 'pc', default = 30, type = 'integer', help = 'number of PC to use [default %(default)s]')
parser$add_argument('--resolution', dest = 'resolution', default = 0.8, type = 'double', help = 'cluster resolution[default %(default)s]')
parser$add_argument('--topN', dest = 'topN', default = 1, type = 'integer', help = 'Top N marker genes in each cluster [default %(default)s]' )

parser$add_argument('--mSpatialVar', dest = 'mSpatialVar', action="store_true", default=FALSE, help = 'Add "--mSpatialVar" to find spatial variable features [default %(default)s]')

opts = parser$parse_args()
print(opts)





###################DEBUG#####################
# ## basic parameters
# rm(list=ls())
# gc()
# opts <- list()
# opts$input <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/qc_test_out/spatial/gem_stereo/ob_out/ob_filt_norm_seuratObject.rds'
# # opts$input <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/qc_test_out/single_cell/10x_brain_out/brain_filt_norm_seuratObject.rds'
# opts$out <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/ob_clsuter1'
# opts$ptsize <- 1
# opts$prefix <- 'brain'
# opts$resolution <- 0.8
# opts$pc = 30
# opts$topN <- 1
#############################################
 


### Function

#' ColorPalette
ColorPalette <- function(number){
    if (number <= 25){
        colorScheme = c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
                        'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
                        'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
                        'green1', 'yellow4', 'yellow3','darkorange4', 'brown')
    } else {
        suppressMessages(library(RColorBrewer))
        qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
        colorScheme <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
        return (colorScheme)
    }
}
#' StereoDimPlot()
StereoDimPlot <- function(object, cols = NULL, pt.size = 1, title = 'Clusters'){
    df <- object@meta.data
    p <- ggplot(df, aes(x = x, y = y, color = seurat_clusters)) +
            geom_point(shape = 19, size = pt.size) +
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
                                  axis.title = element_blank(), axis.line = element_blank(), legend.position = 'right') +
            coord_fixed()  +
            guides(colour = guide_legend(override.aes = list(size=3), nrow = 15, title = title))  +
                  theme_void()

    if (!is.null(cols)){
    p <-  p + scale_color_manual(values = cols)
        }
    return(p)
}

#' Clustering
Clustering <- function(object, dim = 30, resolution = 0.8){
    object <- suppressMessages(RunPCA(object, features = VariableFeatures(object = object)))
    object <- FindNeighbors(object, dims = 1:dim)
    object <- FindClusters(object, resolution = resolution)
    object <- suppressMessages(RunUMAP(object, dims = 1:dim))
    ## check_duplicates这个参数在github的issues里提到，大概是忽略seurat中出现的重复基因/细胞？
    object <- RunTSNE(object, dims = 1:dim, check_duplicates = FALSE)
    return (object)
}

#' Seurat process
DoSeuratProcess <- function(opts){
    dir.create(file.path(opts$out), showWarnings = FALSE)
    setwd(file.path(opts$out))
    
    #' 01 load seurat object
    sce <- readRDS(opts$input)
    # check the variable size
    format(object.size(sce), units = 'auto')

    if (DefaultAssay(sce) == 'SCT'|DefaultAssay(sce) == 'integrated'){
        #' 02 Clustering
        sce <- Clustering(sce, dim = opts$pc, resolution = opts$resolution)
        
        # 03 save intermediate variable
        #' Extract the clusterInfo file
        cluster.info <- as.data.frame(Idents(object = sce))
        colnames(cluster.info) = "cluster"
        write.table(cluster.info, file = paste(opts$out,'/', opts$prefix, '_clusterInfo.txt', sep = ''))
        # save the clusterInfo for cell cell interaction analysis
        cluster.info.cci <- data.frame('cells' = rownames(cluster.info), 'cluster_symbol' =   paste('cluster_', cluster.info$cluster, sep = ''), 'cluster_num' = cluster.info$cluster)
        write.table(cluster.info.cci, file = paste(opts$out,'/', opts$prefix, '_clusterInfo_CCI.txt', sep = ''), sep = '\t', quote=FALSE, row.names = FALSE )

        # save the clusterInfo for TF analysis
        try({
            cluster.info.tf <- data.frame('cells' = rownames(cluster.info), 'cluster' =   paste('cluster_', cluster.info$cluster, sep = ''))
            write.table(cluster.info.tf, file = paste(opts$out,'/', opts$prefix, '_clusterInfo_TF.txt', sep = ''), sep = '\t', quote=FALSE, row.names = FALSE )
        }) 
    
        #' 04 Visualization 
        colorScheme <- ColorPalette(length(levels(sce$seurat_clusters)))
      
        #' 1. visualizing both cells and features that define the PCA 
        pdf(file = paste(opts$out,'/', opts$prefix, '_FeatureDefinePCA.pdf', sep = ''), w = 12, h = 5)
        p <- VizDimLoadings(sce, dims = 1:2, reduction = "pca")
        print (p)
        while (!is.null(dev.list()))  dev.off() 
        
        #' 2. top N PCA DimHeatmap
        pdf(file = paste(opts$out,'/', opts$prefix, '_top15_PCA_DimHeatmap.pdf', sep = ""), w = 10, h = 10)
        p <- DimHeatmap(sce, dims = 1:15, cells = 500, balanced = TRUE)
        print (p)
        while (!is.null(dev.list()))  dev.off()
        
        #' 3. save the PCA, UMAP and TSNE results
        pdf(file = paste(opts$out,'/', opts$prefix, "_PCA.pdf", sep = ""), w = 10, h = 10)
        p <- DimPlot(sce, reduction = "pca", combine = TRUE, cols = colorScheme)
        print (p)
        while (!is.null(dev.list()))  dev.off() 
        
        pdf(file = paste(opts$out,'/', opts$prefix, "_UMAP.pdf", sep = ""), w = 10, h = 10)
        # p <- DimPlot(sce, reduction = "umap",  cols = colorScheme)
        p <-Seurat::DimPlot(sce, reduction = "umap", pt.size = 1.2, label = TRUE, cols = colorScheme, label.size = 8) + 
            theme(text = element_text(size = 20)) + 
            theme(legend.key.height = unit(1.3, 'cm'), legend.text = element_text(size = 25))
        print (p)
        while (!is.null(dev.list()))  dev.off() 
        
        pdf(file = paste(opts$out,'/', opts$prefix, "_TSNE.pdf", sep = ""), w = 10, h = 10)
        p <- DimPlot(sce, reduction = "tsne",   cols = colorScheme)
        print (p)
        while (!is.null(dev.list()))  dev.off() 

        # pdf(file = paste(opts$out,'/', opts$prefix,"_orgin_UMAP.pdf", sep = "")) 
        # DimPlot(sce, reduction = "umap",  group.by  = 'orig.ident') 
        # while (!is.null(dev.list()))  dev.off()  
      
        try({
        
            #' 5. perform cell cycle analysis
            # cc.genes: cell cycle genes
            sce <- CellCycleScoring(sce, s.features = cc.genes$s.genes, g2m.features =cc.genes$g2m.genes, set.ident = F)

            ## Errors occur when visualzing data without cell cycle
            pdf(file = paste(opts$out,'/', opts$prefix, "_cellCycle_UMAP.pdf", sep = ""),w = 10,h = 5)
            # p <- DimPlot(sce, reduction = "umap", group.by  ='Phase', cols=c("yellow", "#458B74", "#483D8B"))
            # print (p)
            p1 <- DimPlot(sce, reduction = "umap",  cols = colorScheme)
            p2 <- DimPlot(sce, reduction = "umap",  group.by  ='Phase', cols=c("yellow", "#458B74", "#483D8B"))
            plot(p1 + p2)
            while (!is.null(dev.list()))  dev.off() 
            
            pdf(file = paste(opts$out,'/', opts$prefix,"_cellCycle_TSNE.pdf", sep = ""), w = 10, h = 5)
            p1 <- DimPlot(sce, reduction = "tsne",  cols = colorScheme)
            p2 <- DimPlot(sce, reduction = "tsne",  group.by  = 'Phase', cols=c("yellow", "#458B74", "#483D8B"))
            plot(p1 + p2)  
            while (!is.null(dev.list()))  dev.off() 
        } )

        #' 6. find markers of each cluster
        # sce.markers <- FindAllMarkers(sce, assay = 'SCT', only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
        if ('RNA' %in% names(sce@assays)){
            assay_used <- 'RNA'
        } else if ('Spatial' %in% names(sce@assays)){
            assay_used <- 'Spatial'
        }
        try({
            DefaultAssay(sce) <- assay_used
            sce <- NormalizeData(sce) 
            FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
            all.genes <- rownames(sce)
            sce <- ScaleData(sce, features = all.genes)

            sce.markers <- FindAllMarkers(sce, assay = assay_used, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25) ## 官方是建议用'RNA'即原始数据进行差异分析。
            # select the marker genes by AUC value
            suppressMessages(library(presto))
            sce.markers.auc <- wilcoxauc(sce, seurat_assay  = assay_used, 'seurat_clusters')
            
            write.table(sce.markers, file = paste(opts$out,'/', opts$prefix, "_markerGene.txt", sep = ""), sep = "\t")
            write.table(sce.markers.auc, file = paste(opts$out,'/', opts$prefix, "_markerGene_auc.txt", sep = ""), sep = "\t")
        })
       
      
        #' 7. visualize the marker gene
        suppressMessages(library(dplyr))
        # Select top N markers of each cluster based on avg_log2FC
        try({
            all.cluster.topn <- sce.markers %>%  group_by(cluster) %>% top_n(n = opts$topN, wt = avg_log2FC)
            write.table(all.cluster.topn, file = paste(opts$out,'/', opts$prefix,"_markerGene_top",opts$topN, ".txt", sep = ""), sep = "\t")
            
            target_gene <- unique(all.cluster.topn$gene)
            if (length(target_gene) > 40){
                target_gene <- target_gene[1:40]
                print("Warning! Only 40 target genes were selected!")
            }
            
            pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_DoHeatmap.pdf', sep = ""),w = 10,h = 10*length(target_gene)/20)
            p <- DoHeatmap(sce, features = target_gene, label = T) +  NoLegend()
            print (p)
            while (!is.null(dev.list()))  dev.off() 
            
            pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_dotPlot.pdf', sep = ""))
            p <- DotPlot(sce, features = target_gene) + RotatedAxis() 
            print (p)
            while (!is.null(dev.list()))  dev.off() 
            
            pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_VinPlo.pdf', sep = ""))
            p <- VlnPlot(sce, features = target_gene, pt.size = 0.5, stack = TRUE) 
            print (p)
            while (!is.null(dev.list()))  dev.off() 
            
            pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_UMAP.pdf', sep = ""),w = 9, h = 4*(length(target_gene)%/%3))
            p <- FeaturePlot(sce, features = target_gene, reduction = "umap", pt.size = 0.5) 
            print (p)
            while (!is.null(dev.list()))  dev.off() 
      
            pdf(file = paste(opts$out,'/', opts$prefix, '_cluster_markers_top',opts$topN,'_FeaturePlot_TSNE.pdf', sep = ""),w = 9, h = 4*(length(target_gene)%/%3))
            p <- FeaturePlot(sce, features = target_gene, reduction = "tsne", pt.size = 0.5) 
            print (p)
            while (!is.null(dev.list()))  dev.off() 
        })
        
        #'Check if images information was included in seurat object (1: spatial; 0: single cell) .
        if (length(sce@images) == 1){
            #' 4. preview the data in Spatial and seurat clusters
            # 图片的排版有问题，暂时分开输出。Seurat绘制空间散点图占内存，需要及时释放??
            pdf(file = paste(opts$out,'/', opts$prefix,"_seurat_clusters_UMAP.pdf", sep = ""), w = 10, h = 5)
            # p1  <- DimPlot(sce, reduction = "umap", cols = colorScheme)
            # p2 <- SpatialPlot(sce, group.by = "seurat_clusters", cols = colorScheme, pt.size.factor = opts$ptsize, stroke = 0)  + theme(legend.position = "right") 
            p2 <- StereoDimPlot(sce, cols = colorScheme, pt.size = opts$ptsize, title = 'seurat_clusters') + theme(legend.position = "right") 
            print(p2)
            while (!is.null(dev.list()))  dev.off() 
            gc()
           
            if (opts$mSpatialVar){
                #' 8. Find Spatial cluster in spatial
                sce <- FindSpatiallyVariableFeatures(sce, assay = "SCT", selection.method = "markvariogram")
                # The cluster info should be included in the following table (Future work)
                sce.markers.st <- SpatiallyVariableFeatures(sce, selection.method = "markvariogram")
                write.table(data.frame('spatialVariableFeature' = sce.markers.st), file = paste(opts$out,'/', opts$prefix, "_spatialVariableFeature.txt", sep = ""), sep = "\t")
                # Draw scatter plots of spatial features separately
                top.features <- head(SpatiallyVariableFeatures(sce, selection.method = "markvariogram"), 6)
                for (i in top.features){
                    pdf(file = paste(opts$out,'/', opts$prefix,"_", i, "_SpatialFeaturesPlots.pdf", sep = ""), w = 5, h = 10)
                    p <- SpatialFeaturePlot(sce, features = i, pt.size.factor = opts$ptsize, alpha = c(0.1, 1), stroke = 0)
                    print (p)
                    while (!is.null(dev.list()))  dev.off() 
                    # release the memory
                    gc()
                }
            }
        }
        
        #' 9. save the data for further analysis
        saveRDS(sce, paste(opts$out, '/', opts$prefix, "_filt_norm_cluster_seuratObject.rds", sep = ""))
    }else {
      stop("The loaded Seurat object should contain SCT assasy", call.=FALSE)
      }
}

### Main Work flow
if ((! is.null(opts$input)) & (! is.null(opts$out)) ){
    # library
    suppressMessages(library(Seurat))
    suppressMessages(library(ggplot2))
    
    DoSeuratProcess(opts)
} else {
    print(parser$print_help())
    stop("You need to provide the following parameters: --input, --out", call.=FALSE)
}
