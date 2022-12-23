args<-commandArgs(T)                                                                                                                                       

# load the parameter
wkd <- args[1]
norm.data <- args[2]
project_name <- args[3]
cluster_resolution <- as.numeric(args[4])
pca.dim.select <- as.numeric(args[5])


#wkd <- "/home/tonyleao/wkd/Scat/tests/liver_st_tests/output/process/cell_cluster"
#norm.data <- "/home/tonyleao/wkd/Scat/tests/liver_st_tests/output/process/normalization/CN73_norm_CCA_C1_trans_filtered_norm.txt.gz"
#project_name <- "CN73_norm_CCA_C1_trans_filtered_norm"
#cluster_resolution <- 1.2 
#pca.dim.select <- 20


# set the work_directory
setwd(wkd)

# Load the data
Ndata <-
  read.table(
    gzfile(norm.data),
    header = T,
    row.name = 1,
    check.names = FALSE
  )

# Configure other parameter
colorScheme = c("#999999", "#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
          "#0072B2", "#D55E00", "#CC79A7", "#990000", "#9900cc", "#66FF66",
          "#663300", "#0000FF", "#CC0033", "#FF0000", "#000099", "#660066",
          "#333333", "#FFCCCC", "#993366", "#33CC33", "#000099", "#CC9900"
)

# output the colorScheme
png(filename = paste(project_name, "colorScheme.png", sep = "_"))
plot(NULL, xlim=c(0,length(colorScheme)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(colorScheme)-1), 0, 1:length(colorScheme), 1, col=colorScheme)
while (!is.null(dev.list()))  dev.off()


# draw histogram
spot_gene=c()
for(i in 2:dim(Ndata)[2]){
  x=Ndata[,i]
  p=which(x>0)
  spot_gene[i-1]=length(p)
}



pdf(
  file = paste(project_name, "gene_dis_per_spot.pdf", sep = "_"),
  w = 4,
  h = 4
)

h <- hist(
  spot_gene,
  breaks = 40,
  xlab = "Number of unique genes per spot",
  ylab = "Number of spots",
  main = ""
)
my <- max(h$counts)
lines(
  x = c(mean(spot_gene), mean(spot_gene)),
  y = c(0, my),
  col = "red",
  lty = 2,
  lwd = 2
)
text(
  x = mean(spot_gene),
  y = my,
  labels = round(mean(spot_gene)),
  col = "red",
  pos = 2
)
while (!is.null(dev.list()))  dev.off()





# Ref to Seurat - Guided Clustering Tutorial
# Link: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
suppressMessages(library(Seurat))
packageVersion("Seurat")
# Initialize the Seurat object with the raw (non-normalized data).
# Create Seurat Object
Ndata.proj <-
  CreateSeuratObject(
    counts = Ndata,
    project = substr(project_name, 1,5),
    min.cells = 3,
    min.features = 200
  )

# QC and selecting cells for further analysis (Not run)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
# I am confused about the pattern in the PercentageFeatureSet(), so I just set pattern = "^MT-"
Ndata.proj[["percent.mt"]] <- PercentageFeatureSet(Ndata.proj, pattern = "^MT-")

p <- suppressWarnings(VlnPlot(Ndata.proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
tiff(filename = paste(project_name, "QC_metrics_violinPlot.tiff", sep = "_"),res = 300,height = 2000,width = 3000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "QC_metrics_violinPlot.pdf", sep = "_"),w = 8,h = 5)
p
while (!is.null(dev.list()))  dev.off() 


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- suppressWarnings(FeatureScatter(Ndata.proj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
plot2 <- FeatureScatter(Ndata.proj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

tiff(filename = paste(project_name, "FeatureScatter.tiff", sep = "_"),res = 300,height = 2000,width = 3000)
plot1 + plot2
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "FeatureScatter.pdf", sep = "_"),w = 8,h = 5)
plot1 + plot2
while (!is.null(dev.list()))  dev.off() 


# Filteration based on the nFeature_RNA and percent.mt (Not run)
# Ndata.proj <- subset(Ndata.proj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)



# Normalizing the data (Not run)
# Ndata.proj <-
#   NormalizeData(Ndata.proj,
#                 normalization.method = "LogNormalize",
#                 scale.factor = 10000)



# Identification of highly variable features (feature selection)
Ndata.proj <-
  FindVariableFeatures(
    Ndata.proj,
    assay = NULL,
    selection.method = "vst",
    nfeatures = 2000
  )


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(Ndata.proj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(Ndata.proj)
plot2 <- suppressMessages(LabelPoints(plot = plot1, points = top10, repel = TRUE))

tiff(filename = paste(project_name, "variableFeaturePlot.tiff", sep = "_"),res = 300,height = 2000,width = 4000)
plot1 + plot2
while (!is.null(dev.list()))  dev.off() 


pdf(file = paste(project_name, "variableFeaturePlot.pdf", sep = "_"),w = 12,h = 5)
plot1 + plot2
while (!is.null(dev.list()))  dev.off() 


# Scaling the data
all.genes <- rownames(Ndata.proj)
Ndata.proj <- ScaleData(Ndata.proj, features = all.genes)

# Perform linear dimensional reduction
Ndata.proj <-
  suppressMessages(RunPCA(Ndata.proj, features = VariableFeatures(object = Ndata.proj)))

# Visualizing both cells and features that define the PCA 
p <- VizDimLoadings(Ndata.proj, dims = 1:2, reduction = "pca")
tiff(filename = paste(project_name, "featureDefinePCA.tiff", sep = "_"),res = 300,height = 2000,width = 4000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "featureDefinePCA.pdf", sep = "_"),w = 12,h = 5)
p
while (!is.null(dev.list()))  dev.off() 

# pca
p <- DimPlot(Ndata.proj, reduction = "pca", split.by = 'ident')
tiff(filename = paste(project_name, "PCA.tiff", sep = "_"),res = 300,height = 2000,width = 2000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "PCA.pdf", sep = "_"),w = 10,h = 10)
p
while (!is.null(dev.list()))  dev.off() 

# top15_DimHeatmap
tiff(filename = paste(project_name, "top15_DimHeatmap.tiff", sep = "_"),res = 300,height = 4000,width = 4000)
DimHeatmap(Ndata.proj, dims = 1:15, cells = 500, balanced = TRUE)
while (!is.null(dev.list()))  dev.off()

pdf(file = paste(project_name, "top15_DimHeatmap.pdf", sep = "_"),w = 10,h = 10)
DimHeatmap(Ndata.proj, dims = 1:15, cells = 500, balanced = TRUE)
while (!is.null(dev.list()))  dev.off()



# Determine the 'dimensionality' of the dataset
# NOTE: This process can take a long time for big datasets, comment out for expediency. 
# More approximate techniques such as those implemented in ElbowPlot() can be used to 
# reduce computation time.
dim_max <- length(Ndata.proj@reductions$pca)

if (pca.dim.select > dim_max) {
  pca.dim.select <- dim_max
  print ("Warning! The pca.dim.select exceeds the pca dimention calculated!")
  sprintf('pca.dim.select has been set to %s',pca.dim.select)
}


Ndata.proj <- JackStraw(Ndata.proj, num.replicate = 100)
Ndata.proj <- ScoreJackStraw(Ndata.proj, dims = 1:20)

# select significant PCA dimension (p.values < 0.05) (Not Run)
# Here we configure the parameter "pca.dim.select" by the output
# overall.p.values <- Ndata.proj@reductions$pca@jackstraw$overall.p.values
# overall.p.values <- as.data.frame(overall.p.values)
# select.p.values <- overall.p.values[overall.p.values$Score < 0.05,]
# pca.dim.select <- nrow(select.p.values)
# pca.dim.select <- select.p.values$PC[nrow(select.p.values)]


# if (nrow(overall.p.values)==nrow(select.p.values)) {
#   print("Warning! The first 20 PCA components may not be sufficient to represent all the data.")
# }


# Output the JackStrawPlot and ElbowPlot
p1 <- JackStrawPlot(Ndata.proj, dims = 1:20)
p2 <- ElbowPlot(Ndata.proj)

tiff(filename = paste(project_name, "JackStrawPlot_ElbowPlot.tiff", sep = "_"),res = 300,height = 2000,width = 4000)
p1 + p2
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "JackStrawPlot_ElbowPlot.pdf", sep = "_"),w = 10,h = 5)
p1 + p2
while (!is.null(dev.list()))  dev.off()



# Cluster the cells/spot
pca.dim.select
Ndata.proj  <- FindNeighbors(Ndata.proj, dims = 1:pca.dim.select)
Ndata.proj <- FindClusters(Ndata.proj, resolution = cluster_resolution)
# head(Idents(Ndata.proj), 5)
cluster.info <- Idents(object = Ndata.proj)
# head(cluster.info)


write.table(
  cluster.info,
  file = paste(project_name, "clusterInfo.txt", sep = "_"),
  sep = "\t",
  quote = FALSE
)


# Run non-linear dimensional reduction (UMAP/tSNE)
Ndata.proj <- RunUMAP(Ndata.proj, dims = 1:pca.dim.select)
Ndata.proj <- RunTSNE(Ndata.proj, dims = 1:pca.dim.select)


p <- DimPlot(Ndata.proj, reduction = "umap",
             # group.by = "orig.ident",
             # group.by = "nFeature_RNA",
             cols = colorScheme[1:length(levels(Ndata.proj$seurat_clusters))])
tiff(filename = paste(project_name, "UMAP.tiff", sep = "_"),res = 300,height = 2000,width = 2000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "UMAP.pdf", sep = "_"),w = 10,h = 10)
p
while (!is.null(dev.list()))  dev.off() 


p <- DimPlot(Ndata.proj, reduction = "tsne",
             # group.by = "orig.ident",
             # group.by = "nFeature_RNA",
             cols = colorScheme[1:length(levels(Ndata.proj$seurat_clusters))])
tiff(filename = paste(project_name, "TSNE.tiff", sep = "_"),res = 300,height = 2000,width = 2000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "TSNE.pdf", sep = "_"),w = 10,h = 10)
p
while (!is.null(dev.list()))  dev.off() 


# Finding differentially expressed features (cluster biomarkers)
# FindMarkers() for makers in specific cluster, e.g. 1, 2, 3
# Usage:
# cluster2.markers <- FindMarkers(Ndata.proj, ident.1 = 2, min.pct = 0.25)
# head(cluster2.markers, n = 5)


Ndata.proj.markers <-
  FindAllMarkers(
    Ndata.proj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )

write.table(
  Ndata.proj.markers,
  file = paste(project_name, "markerGene.txt", sep = "_"),
  sep = "\t"
)
            

# Visualize the marker gene
suppressMessages(library(dplyr))
# Select top 5 markers of each cluster based on avg_log2FC
all.cluster.top5 <-
  Ndata.proj.markers %>%  group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

write.table(
  all.cluster.top5,
  file = paste(project_name, "top5_markerGene.txt", sep = "_"),
  sep = "\t"
)


target_gene <- unique(all.cluster.top5$gene)
if (length(target_gene) > 40){
  target_gene <- target_gene[1:40]
  print("Warning! Only top 40 target genes were selected!")
}


p <- DoHeatmap(Ndata.proj,
               features = target_gene,
               label = T) +  NoLegend()
tiff(filename = paste(project_name, "top5_markers_DoHeatmap.tiff", sep = "_"),res = 300,height = 2000*length(target_gene)/20,width = 2000)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "top5_markers_DoHeatmap.pdf", sep = "_"),w = 10,h = 10*length(target_gene)/20)
p
while (!is.null(dev.list()))  dev.off() 


p <- DotPlot(Ndata.proj, features = target_gene) + RotatedAxis() 

tiff(filename = paste(project_name, "top5_markers_dotPlot.tiff", sep = "_"),res = 300,height = length(levels(Ndata.proj))*200,width = length(target_gene)*110)
p
while (!is.null(dev.list()))  dev.off() 

pdf(file = paste(project_name, "top5_markers_dotPlot.pdf", sep = "_"),w = 0.4*length(target_gene),h = 0.65*length(levels(Ndata.proj)))
p
while (!is.null(dev.list()))  dev.off() 


p <- FeaturePlot(
  Ndata.proj,
  features = target_gene ,
  reduction = "umap",
  pt.size = 0.5
) 


tiff(filename = paste(project_name, "top5_markers_FeaturePlot_UMAP.tiff", sep = "_"),res = 300,height = length(target_gene)%/%4*1000,width = 4000)
p
while (!is.null(dev.list()))  dev.off() 
pdf(file = paste(project_name, "top5_markers_FeaturePlot_UMAP.pdf", sep = "_"),w = 10,h = 2.5*(length(target_gene)%/%4))
p
while (!is.null(dev.list()))  dev.off() 



p <- FeaturePlot(
  Ndata.proj,
  features = target_gene ,
  reduction = "tsne",
  pt.size = 0.5
) 

tiff(filename = paste(project_name, "top5_markers_FeaturePlot_TSNE.tiff", sep = "_"),res = 300,height = length(target_gene)%/%4*1000,width = 4000)
p
while (!is.null(dev.list()))  dev.off() 
pdf(file = paste(project_name, "top5_markers_FeaturePlot_TSNE.pdf", sep = "_"),w = 12,h = 2.5*(length(target_gene)%/%4))
p
while (!is.null(dev.list()))  dev.off() 


p <- VlnPlot(
  Ndata.proj,
  features = target_gene,
  pt.size = 0.5,
  stack = TRUE,
) 

tiff(filename = paste(project_name, "top5_markers_VinPlo.tiff", sep = "_"),res = 300,height = length(levels(Ndata.proj))*500,width = length(target_gene)%/%4*1000)
p
while (!is.null(dev.list()))  dev.off()
pdf(file = paste(project_name, "top5_markers_VinPlo.pdf", sep = "_"),w = 0.6*length(target_gene),h = 0.2*length(target_gene))
p
while (!is.null(dev.list()))  dev.off() 


# Assigning cell type identity to clusters (Not ru n)

# new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK")
# names(new.cluster.ids) <- levels(Ndata.proj)
# Ndata.proj <- RenameIdents(Ndata.proj, new.cluster.ids)
# 
# 
# p <- DimPlot(Ndata.proj, reduction = "umap", 
#         cols = colorScheme[1:length(levels(Ndata.proj$seurat_clusters))],
#         label = TRUE, pt.size = 0.5) 
# 
# tiff(
#   filename = paste(project_name, "cellTye_UMAP.tiff", sep = "_"),
#   res = 300,
#   height = 2000,
#   width = 2500
# )
# p
# dev.off()



# Save the data for further analysis
saveRDS(Ndata.proj,
        paste(project_name, "seuratObject.rds", sep = "_"))
