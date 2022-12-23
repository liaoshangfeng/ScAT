##单细胞rds转表达矩阵
inputfile='/ldfssz1/ST_OCEAN/USER/liaoshangfeng/Project/01.scat_debug/data/mouse.brain/Sc.mouse.brain.rds'
outfile='sc.mouse.brain.tsv'

data=readRDS(inputfile)
df=data.frame(data@assays$RNA@counts)
write.table(df,file=outfile,quote=F,sep='\t')



