###空间rds转表达矩阵
##修改输入和输出文件名称

inputfile='/ldfssz1/ST_OCEAN/USER/liaoshangfeng/Project/01.scat_debug/data/mouse.brain/ST.mouse.brain.rds'
outfile='./ST.ms.brain.tsv'



Stdata=readRDS(inputfile,assay = "Spatial", verbose = FALSE)
coord=GetTissueCoordinates(
  Stdata,
  scale = "lowres",
  cols = c("imagerow", "imagecol"),
)
coord$imagerow<-as.character(coord$imagerow)
coord$imagecol<-as.character(coord$imagecol)
df_coord=tidyr::unite(coord, "imagerow_imagecol", imagerow, imagecol)
expr_raw <- GetAssayData(object = anterior, assay.type = "Spatial", slot = "counts")
expr <- as(Class = 'matrix', object = expr_raw)
colnames(expr) <- df_coord[[1]]
write.table(expr,file=outfile,quote=F,sep='\t')


