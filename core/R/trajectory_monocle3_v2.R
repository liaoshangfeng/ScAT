## Trajectory analysis by monocle3

parser = argparse::ArgumentParser(description = 'Script for trajectory analysis')
parser$add_argument('-i', '--input', dest = 'input', default=NULL, help = 'SeuratObject RDS file [default %(default)s]')
parser$add_argument('-o', '--out', dest = 'out', default=NULL, help = 'Directory to save file [default %(default)s]')
parser$add_argument('-g', '--gene', dest = 'gene', default=NULL, help = 'Gene name to define root node of pseudotime trajectory [default %(default)s]')
parser$add_argument('--mDEG', dest = 'mDEG', action='store_true', default=FALSE, help = 'Allow differential expression analysis [default %(default)s]')
parser$add_argument('--mGCN', dest = 'mGCN', action="store_true", default=FALSE, help = 'Allow gene co-expression analysis [default %(default)s]')
parser$add_argument('--mTop_markers', dest = 'mTop_markers', action="store_true", default=FALSE, help = 'Allow top markers analysis [default %(default)s]')
parser$add_argument('--prefix', dest = 'prefix', default = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident [default %(default)s]')

opts = parser$parse_args()
# print(opts)

###################DEBUG#####################
# rm(list=ls())
# gc()
# 
# opts <- list()
# # opts$input <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/brain_clsuter1/brain_filt_norm_cluster_seuratObject.rds'
# # opts$out <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/traj_brain'
# opts$input <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/ob_clsuter1/ob_filt_norm_cluster_seuratObject.rds'
# opts$out <- '/home/tonyleao/wkd/dev/scat_dev/cluster_module/traj_ob1/'
# # opts$input <- "/home/tonyleao/wkd/data_to_viz/scat_dev/cluster_module/cluster_out/ob_filt_norm_cluster_seuratObject.rds"
# # opts$out <- "/home/tonyleao/wkd/data_to_viz/scat_dev/cluster_module/traj_out"
# opts$gene <- NULL
# opts$mTop_markers <- TRUE
# opts$mGCN <- TRUE
# opts$mDEG <- TRUE
# opts$prefix <- 'ob'
#############################################


### Function

#' define function to select root nodes
get_earliest_principal_node <- function(cds, time_bin = '0'){
    #  select the cluster 0 as root node
    cell_ids <- which(colData(cds)[, 'seurat_clusters'] == time_bin)
    closest_vertex <- cds@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[['UMAP']])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
}

#' define main trajectory process 
trajectory_process <- function(opts){
    # 01 Import the seurat object rds file
    sce <-  readRDS(opts$input)
    # names(so@assays) = c("RNA", "SCT")
    
    cds <- as.cell_data_set(sce)
    ## Calculate size factors using built-in function in monocle3
    cds <- estimate_size_factors(cds)
    ## Add gene names into CDS
    cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sce[[DefaultAssay(sce)]])

    # if (DefaultAssay(sce) == 'SCT'){
    #     ## for spatial seurat object or seurat object undergo SCT normalization
    #     cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sce[["SCT"]])
    # } else if (DefaultAssay(sce) == 'RNA'){
    #     cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sce[["RNA"]])
    # } else if (DefaultAssay(sce) == 'Spatial'){
    #     cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sce[["Spatial"]])
    # } else if (DefaultAssay(sce) == 'integrated'){
    #     cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(sce[["integrated"]])
    # }



    
    ## Set reduction method
    cds <- cluster_cells(cds,reduction_method = "UMAP")
    # 02 Identify the trajectory
    cds <- learn_graph(cds)
    # 03 order cells
    if (! is.null(opts$gene)){
        print(opts$gene)
        root_pr_nodes = opts$gene
    } else {
        root_pr_nodes = get_earliest_principal_node(cds)
        print(opts$gene)
    }
    cds <- order_cells(cds, root_pr_nodes=root_pr_nodes)
    return(list(sce, cds))
}

#' define function for differential expression analysis
diff_express_analysis <- function(cds, opts){
    diff_gene <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
    id <- row.names(subset(diff_gene, q_value < 0.05))
    write.csv(id, paste(opts$out, '/', opts$prefix, "_trajectory_diff_genes.csv",sep = ""), row.names = F)
    write.csv(diff_gene, paste(opts$out, '/', opts$prefix,  "_trajectory_all_genes.csv",sep = ""), row.names = T)
    
    diff_gene$gene_short_name <- rownames(diff_gene)
    track_genes_sig <- diff_gene %>% 
      filter(q_value < 0.05) %>% 
      top_n(n=10, morans_I) %>%
      pull(gene_short_name) %>% 
      as.character()
    
    p3 <- plot_genes_in_pseudotime(cds[track_genes_sig,], color_cells_by="pseudotime", min_expr=0.5, ncol = 2)
    ggsave(paste(opts$out, '/', opts$prefix, "_genes_Jitterplot_pseudotime.pdf",sep=""), plot = p3, width = 8, height = 6)
    
    p4 <- plot_genes_in_pseudotime(cds[track_genes_sig,], color_cells_by="seurat_clusters", min_expr=0.5, ncol = 2)
    ggsave(paste(opts$out, '/', opts$prefix, "_genes_Jitterplot_cluster.pdf",sep=""), plot = p4, width = 8, height = 6)
    # Feature plot
    p5 <- plot_cells(cds, genes=track_genes_sig, show_trajectory_graph=FALSE,
                     label_cell_groups=FALSE,  label_leaves=FALSE)
    p5$facet$params$ncol <- 5
    ggsave(paste(opts$out, '/', opts$prefix, "_genes_Featureplot.pdf", sep = ""), plot = p5, width = 20, height = 8)
    return(diff_gene)
}

#' define function for exploring co-expression module
find_co_expression_module <- function(diff_gene, cds,  opts){
    genelist <- pull(diff_gene, gene_short_name) %>% as.character()
    gene_module <- find_gene_modules(cds[genelist,], resolution=1e-1, cores = 6)
    write.csv(gene_module, paste(opts$out, '/', opts$prefix, "_genes_module.csv",sep=""), row.names = F)
    cell_group <- tibble::tibble(cell=row.names(colData(cds)), 
                                 cell_group=colData(cds)$seurat_clusters)
    agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
    row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
    p6 <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
    ggsave(paste(opts$out, '/', opts$prefix, "_genes_module.pdf",sep=""), plot = p6, width = 8, height = 8)
    saveRDS(agg_mat, file = paste(opts$out, '/', opts$prefix, "_aggregate_gene_expression.rds", sep="")) 
}

#' define function for identifing top markers
find_top_markers <- function(cds, opts){
    marker_test_res = top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=8)
    write.csv(marker_test_res, paste(opts$out, '/', opts$prefix, "_top_markers.csv",sep=""), row.names = F)
    
    top_specific_markers = marker_test_res %>%
      filter(fraction_expressing >= 0.10) %>%
      group_by(cell_group) %>%
      top_n(1, pseudo_R2)
    
    top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
    
    p8<- plot_genes_by_group(cds,
                        top_specific_marker_ids,
                        group_cells_by="partition",
                        ordering_type="maximal_on_diag",
                        max.size=3)
    ggsave(paste(opts$out, '/', opts$prefix, "_top1_marker_gene_partition.pdf",sep=""), plot = p8, width = 8, height = 8)
    
    
    top_specific_markers = marker_test_res %>%
      filter(fraction_expressing >= 0.10) %>%
      group_by(cell_group) %>%
      top_n(3, pseudo_R2)
    
    top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
    
    p9 <- plot_genes_by_group(cds,
                        top_specific_marker_ids,
                        group_cells_by="partition",
                        ordering_type="cluster_row_col",
                        max.size=3)
    ggsave(paste(opts$out, '/', opts$prefix, "_top3_marker_gene_partition.pdf",sep=""), plot = p9, width = 8, height = 8)
    
    
    p10 <- plot_genes_by_group(cds,
                        top_specific_marker_ids,
                        group_cells_by="seurat_clusters",
                        ordering_type="cluster_row_col",
                        max.size=3)
    ggsave(paste(opts$out, '/', opts$prefix, "_top3_marker_gene_seurat_clusters.pdf",sep=""), plot = p10, width = 8, height = 8)
}

main <- function(opts){
    dir.create(file.path(opts$out), recursive = TRUE, showWarnings = FALSE)
    setwd(file.path(opts$out))
    
    # 01  trajectory process
    results <- trajectory_process(opts)
    sce <- results[[1]]
    cds <- results[[2]] 
  
    p <- plot_cells(cds, 
               color_cells_by = 'seurat_clusters',
               label_groups_by_cluster = TRUE, 
               label_leaves = FALSE, 
               label_branch_points = TRUE,
               label_cell_groups = FALSE,
               # graph_label_size = 1.5,
               cell_stroke = 0.8)
    ggsave(paste(opts$out, '/', opts$prefix, "_trajectory.pdf", sep=""), plot = p, width = 5, height = 5)

    p1 <- plot_cells(cds, 
               color_cells_by = 'seurat_clusters',
               label_groups_by_cluster = TRUE, 
               label_leaves = FALSE, 
               label_branch_points = TRUE,
               label_cell_groups = FALSE,
               # graph_label_size = 1.5,
               cell_stroke = 0.8)
    ggsave(paste(opts$out, '/', opts$prefix, "_trajectory_branch.pdf", sep=""), plot = p1, width = 5, height = 5)

    p2 <- plot_cells(
                cds,
                color_cells_by = "pseudotime",
                label_cell_groups = FALSE,
                label_leaves = FALSE,
                label_branch_points = TRUE,
                graph_label_size = 1.5,
                cell_stroke = 0.8)
    ggsave(paste(opts$out, '/', opts$prefix, "_trajectory_pseudotime.pdf", sep=""), plot = p2, width = 5, height = 5)
   
    # 02 differential expression analysis and identification of co-expreession module
    if (opts$mDEG){     
        diff_gene <- diff_express_analysis(cds, opts)
        if (opts$mGCN){
            find_co_expression_module(diff_gene, cds, opts)
        }
    }

    # 03 identify top markers
    if (opts$mTop_markers & length(levels(cds@clusters$UMAP$partitions)) >1){
        find_top_markers(cds, opts)
    }

    # 04 save pseudotime results
    # Extract the results of the quasi - time analysis and return the Seurat object
    pseudotime <- pseudotime(cds, reduction_method = 'UMAP')
    pseudotime <- pseudotime[rownames(sce@meta.data)]
    sce$pseudotime <- pseudotime
    saveRDS(sce, file = paste(opts$out, '/', opts$prefix, "_sce_pseudotime.rds",sep=""))

}


### Main Work flow
if ((! is.null(opts$input))  & (!is.null(opts$out))){
    suppressMessages(library(monocle3))
    suppressMessages(library(Seurat))
    suppressMessages(library(SeuratWrappers))
    suppressMessages(library(patchwork))
    suppressMessages(library(ggplot2))
    suppressMessages(library(dplyr))
    main(opts)
} else {
    print(parser$print_help())
    stop("You need to provide the following parameters: --input, --out", call.=FALSE)
}
