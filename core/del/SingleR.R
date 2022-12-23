args<-commandArgs(T)                                                                                                                                       

# load the parameter
wkd <- args[1]
norm.data <- args[2]
ref_db <- args[3] 
ref_db_col <- args[4]
cell_type <- args[5] 
cluster_file <- args[6]

project_name <- unlist(strsplit(basename(norm.data),split=".",fixed=TRUE))[1]

# set the work_directory
setwd(wkd)
getwd()
print (project_name)


# Configure other parameter
colorScheme = c("#999999", "#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#990000", "#9900cc", "#66FF66",
                "#663300", "#0000FF", "#CC0033", "#FF0000", "#000099", "#660066",
                "#333333", "#FFCCCC", "#993366", "#33CC33", "#000099", "#CC9900"
)

# output the colorScheme
tiff(
    filename = paste(project_name, "colorScheme.tiff", sep = "_"),
    res = 300,
    height = 2000,
    width = 3000
)
plot(NULL, xlim=c(0,length(colorScheme)), ylim=c(0,1), 
     xlab="", ylab="", xaxt="n", yaxt="n")
rect(0:(length(colorScheme)-1), 0, 1:length(colorScheme), 1, col=colorScheme)
while (!is.null(dev.list()))  dev.off()




# load expression data
Ndata <-
    read.table(
        gzfile(norm.data),
        header = T,
        row.names = 1 ,
        check.names = FALSE)




## ---- echo=FALSE, results="hide", message=FALSE-------------------------------
knitr::opts_chunk$set(error=FALSE, message=FALSE, warning=FALSE)
suppressMessages(library(BiocStyle))
suppressMessages(library(SingleR))
suppressMessages(library(ExperimentHub))
suppressMessages(library(DelayedArray))


## Function
.rm_NAs <- function(mat, rm.NA = "rows"){
    # Identify them first before removal, to ensure that 
    # the same rows are columns are removed with 'both'.
    if(rm.NA == "rows" || rm.NA == "both"){
        keep_rows <- !rowAnyNAs(DelayedArray(mat))
    } else {
        keep_rows <- !logical(nrow(mat))
    }
        
    if (rm.NA=="cols" || rm.NA == "both") {
        keep_cols <- !colAnyNAs(DelayedArray(mat))
    } else {
        keep_cols <- !logical(ncol(mat))
    }
    # Avoid making unnecessary copies if possible.
    if (!all(keep_rows)) {
        mat <- mat[keep_rows,,drop=FALSE]
    }
    if (!all(keep_cols)) {
        mat <- mat[,keep_cols,drop=FALSE]
    }
    mat
}


# perform the cell type annotation
all.assays <- list()
reference <- readRDS(ref_db)
all.assays[['logcounts']]<-.rm_NAs(reference, c("rows","cols","both","none"))
ref_colname <- readRDS(ref_db_col)

arg <- list()
arg$colData <- ref_colname
ref_se <- do.call(SummarizedExperiment, c(list(assays=all.assays),arg))

pred1 <- SingleR(test = Ndata, ref = ref_se, assay.type.test=1, labels = ref_se$label.main)
pred2 <- SingleR(test = Ndata, ref = ref_se, assay.type.test=1, labels = ref_se$label.fine)
pred<-cbind(rownames(pred1), pred1$labels, pred2$labels)
colnames(pred)<-c('cells','main','fine')
write.table(pred, cell_type, sep='\t', quote=FALSE, row.names = FALSE)   




# visualize the cell component in cluster

# create the cell_type data.frame
cell_type <- as.data.frame(pred)
rownames(cell_type) <- cell_type[,1]
cell_type <- cell_type[,-1]
colnames(cell_type) = c("main", "fine")

# load cell cluster file
cell_cluster <-
    read.table(
        cluster_file,
        row.names = 1,
        header = T,
        check.names = F
    )
## merge the celltype and clusterInfo
anno_cell <- intersect(rownames(cell_type), rownames(cell_cluster))
df <- cell_type[rownames(cell_type) %in% anno_cell, ]

df$cluster <- rep(NA, nrow(df))
for (cell in anno_cell) {
    df[cell, ]$cluster  = cell_cluster[cell, ]
}

# save the cluster cell component file
write.table(
    df,
    paste(project_name, 'cluter_cell_component.txt', sep = '_'),
    sep = '\t',
    quote = FALSE,
    row.names = FALSE
)



suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))

df$number <- 1

tiff(
    filename = paste(project_name,"cluster_cell_component.tiff", sep = "_"),
    width = 400*length(unique(df$cluster)),
    height = 2500,
    res = 300
)
p <- ggplot(df, aes(cluster, number, fill = main)) +
    geom_bar(stat = "identity", position = "stack") +
    # coord_flip() +
    theme_bw() +
    # scale_fill_wsj("colors6", "")+
    scale_fill_manual(values = colorScheme[1:length(unique(df$main))]) +
    theme(axis.ticks.length = unit(0.5, 'cm')) +
    guides(fill = guide_legend(title = NULL))


p <- p + theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
) + labs(x = "Cell Number", y = "Cell Cluter", title = "") + theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(
        hjust = 1,
        vjust = -40,
        size = 14
    )
)
p
while (!is.null(dev.list()))  dev.off()




# Draw a percentage stack diagram
suppressMessages(library(plyr))
df <- ddply(df, "cluster", transform,
      percent_weight = number / sum(number) * 100)

# library(forcats) 
# order the percentage weight demo
# TotalBySource$Source <- fct_rev(TotalBySource$Source) 

tiff(
    filename = paste(project_name,"cluster_cell_component_percentage.tiff", sep = "_"),
    width = 400*length(unique(df$cluster)),
    height = 2500,
    res = 300
)
p <- ggplot(df, aes(cluster, percent_weight, fill = main)) +
    geom_bar(stat = "identity", position = "stack") +
    # coord_flip() +
    theme_bw() +
    # scale_fill_wsj("colors6", "")+
    scale_fill_manual(values = colorScheme[1:length(unique(df$main))]) +
    theme(axis.ticks.length = unit(0.5, 'cm')) +
    guides(fill = guide_legend(title = NULL))


p <- p + theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
) + labs(x = "Cell Cluster", y = "Cell component", title = "") + theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(
        hjust = 1,
        vjust = -40,
        size = 14
    )
)
p
while (!is.null(dev.list()))  dev.off()
