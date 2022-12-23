library(Seurat)
#for(filedir in c('3536-1-210901','3536-2-210901','3532-2-210901','3532-1-210901')) {

 #counts <- Read10X(data.dir = filedir, gene.column = 1)
 #   sample_rds <- CreateSeuratObject(counts,project=filedir,min.cells=3)
 #   saveRDS(sample_rds,file=paste(filedir,'.rds',sep=''))  
 # }


inputfile='/jdfssz1/ST_TSCBI/PROJECT_temp/USER/huangke/20.muscle/muscle_10matrix.txt'
#输入文件三列
#1.project name
#2.样本名称(示例文件为技术重复)
#3.三个测序文件存放的文件夹

   data=read.table(inputfile,header=FALSE)
   for(sampleid in unique(data$V1)) { 
   # sampleid='TQD210612554'
    print(sampleid) 
      scRNAlist = {}
   for (i in (1:length(subset(data,data$V1==sampleid)$V3))){
    print(subset(data,data$V1==sampleid)$V3[i])
    counts <- Read10X(data.dir = subset(data,data$V1==sampleid)$V3[i],  gene.column = 1)
    scRNAlist[[i]] <- CreateSeuratObject(counts,project=sampleid)
    scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id = subset(data,data$V1==sampleid)$V2[i] )
    }
    print(scRNAlist)
    scRNA=scRNAlist[[1]]
    for (j in (2:length(scRNAlist))){
    #scRNA=scRNAlist[[1]]
    scRNA <- merge(scRNA,scRNAlist[j])
    saveRDS(scRNA,file=paste(sampleid,'rds',sep='.'))
    
    }
    }


