args <- commandArgs(T)

# load the parameter
wkd <- args[1]
non.norm.data <- args[2]
norm.data.name <- args[3]
norm_method <- args[4] 


# set the work_directory
setwd(wkd)


## define the normalization function
norm_linnorm <- function(Ndata,  norm.data.name){
  suppressMessages(library(Linnorm))
  norm.data <-  Linnorm.Norm(Ndata, output = "XPM")
  f_out <- gzfile(norm.data.name, "w")
  write.table(norm.data, f_out, sep = "\t")
  close(f_out)
}


norm_scran <- function(Ndata, norm.data.name){
  suppressMessages(library(scran))
  sce <- SingleCellExperiment(list(counts=Ndata))
  sce <- scran::computeSumFactors(sce)
  sce <- scater::logNormCounts(sce)
  f_out <- gzfile(norm.data.name, "w")
  write.table(sce@assays@data@listData$logcounts, f_out, sep = "\t")
  close(f_out)
}


# load the data
Ndata <-
  read.table(
    gzfile(non.norm.data),
    header = T,
    row.names = 1 ,
    check.names = FALSE)

# excute normalization based on norm_method
if (norm_method == "Linnorm"){
  print ("Here we excute the normalization by Linnorm.")
  norm_scran(Ndata,  norm.data.name)
  
}else if (norm_method == "scran"){
  print ("Here we excute the normalization by scran.")
  norm_scran(Ndata,  norm.data.name)
}
