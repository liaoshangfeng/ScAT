### Get the parameters
library(argparse, quietly = TRUE)
parser = argparse::ArgumentParser(description = 'Spatial pattern identification by SPARK')
parser$add_argument('-i', '--input', dest = 'input', help = 'input seurat object *.rds filename')
parser$add_argument('-o', '--out', dest = 'out', help = 'directory where to save the output files, all output files will be indexed by sample ID')
parser$add_argument('--prefix', dest = 'prefix', default = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident [default %(default)s]')
opts = parser$parse_args()
print(opts)



## Functions
run_SPARK <- function(opts){
    sce <- readRDS(opts$input)
    rawcount <- sce@assays$Spatial@counts
    rawcount[1:5,1:5]
    df <- sce@meta.data
    head(df)
    rm(sce)
    gc()


    info <- data.frame('x' = df$x,
                      'y' = df$y,
                      'total_counts' = df$nCount_Spatial)

    head(info)

    rownames(info) <- rownames(df)
    identical(rownames(df), colnames(rawcount))

    ## Create a SPARK object
    tic('Create a SPARK object:')
    ## filter genes and cells/spots and 
    spark <- CreateSPARKObject(counts=rawcount, 
                                 location=info[,1:2],
                                 percentage = 0.1, 
                                 min_total_counts = 10)
    toc()

    ## total counts for each cell/spot
    spark@lib_size <- apply(spark@counts, 2, sum)

    ## Estimating Parameter Under Null
    tic('Estimating Parameter Under Null:')
    spark <- spark.vc(spark, 
                       covariates = NULL, 
                       lib_size = spark@lib_size, 
                       num_core = 5,
                       verbose = F)
    toc()

    ## Calculating pval
    tic('Calculating pval:')
    spark <- spark.test(spark, 
                         check_positive = T, 
                         verbose = F)
    toc()
    saveRDS(spark, paste(opts$out, '/', opts$prefix, '_spark_Obj.rds'))

    results <- spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
    write.csv(results, paste(opts$out, '/', opts$prefix, '_spark_results.csv'), quote =F)
}

### Main Work flow
if ((! is.null(opts$input)) & (! is.null(opts$out)) ){
    # loading library
    suppressMessages({
    library(SPARK)
    library(Seurat)
    library(tictoc)   
    })
    
    run_SPARK(opts)
    
} else {
    print(parser$print_help())
    stop("You need to provide the following parameters: --input, --out", call.=FALSE)
}