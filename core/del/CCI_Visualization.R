args<-commandArgs(T)                                                                                                                                       
# load the parameter
wkd <- args[1]
count_file <- args[2]
p_value_file <- args[3]
project_name <- args[4]

#wkd <- "/home/tonyleao/wkd/viz_script/cell_cell_interaction_circors/"
#count_file <- "counts.txt"
#p_value_file <- "pvalue.txt"
#project_name <- "CN73_norm_CCA_C1_trans"


# set the work_directory
setwd(wkd)
# Configure other parameter
# colorScheme = c("#999999", "#FF0099", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
#                 "#0072B2", "#D55E00", "#CC79A7", "#990000", "#9900cc", "#66FF66",
#                 "#663300", "#0000FF", "#CC0033", "#FF0000", "#000099", "#660066",
#                 "#333333", "#FFCCCC", "#993366", "#33CC33", "#000099", "#CC9900"
# )

colorScheme <- c("#40A4D8", "#33BEB7", "#B2C224", "#FECC2F", "#FBA127", "#F66320",
                             "#DB3937", "#A463D7", "#0C5BCE")
                 

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



# draw chord plot 
# load the cell and cell interaction data (count)
myData <- read.table(count_file, sep = "\t", header = T, check.names = F)

## format the data to three line (from, to, val)
d <- data.frame("from"=character(0), 
                "to"=character(0),
                "val"=numeric(0))

for (i in colnames(myData)){
  d1 <- data.frame("from"=rownames(myData), "to"=rep(i, nrow(myData)), "cor"=myData[,i])
  d <- rbind(d, d1)
}


d2 <- data.frame("from" = unique(d$from),
                 "cols" =  colorScheme[1:length(unique(d$from))],
                 "order" = seq(1:length(unique(d$from))),
                 "reg1" = unique(d$from))


suppressMessages(library(circlize))
suppressMessages(library(dplyr))
tiff(
  filename = paste(project_name, "_CCI_chord.tiff", sep = "_"),
  width = 2500,
  height = 2500,
  res = 300
)

chordDiagram(
  x = d,
  directional = 1,
  order = d2$from,
  grid.col = d2$cols,
  annotationTrack = "grid",
  transparency = 0.25,
  annotationTrackHeight = c(0.05, 0.1),
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  diffHeight  = -0.04
)

circos.track(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    reg1 = d2 %>% filter(from == sector.index) %>% pull(reg1)
    # reg2 = d1 %>% filter(from == sector.index) %>% pull(reg2)
    # circos.text(x = mean(xlim), y = ifelse(is.na(reg2), 3, 4),labels = reg1, facing = "bending", cex =0.8)
    circos.text(
      x = mean(xlim),
      y = 2.75,
      labels = reg1,
      facing = "downward",
      cex = 0.8
    )
    circos.axis(
      h = "top",
      labels.cex = 0.6,
      labels.niceFacing = FALSE,
      labels.pos.adjust = FALSE
    )
  }
)
while (!is.null(dev.list()))  dev.off()



# draw heatmap
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))

heat <- read.table(file = p_value_file, header = T)
heat_matrix <- as.matrix(heat)

tiff(
  filename = paste(project_name, "CCI_heatmap.tiff", sep = "_"),
  width = 2500,
  height = 2500,
  res = 300
)
pheatmap(heat_matrix)
while (!is.null(dev.list()))  dev.off()
