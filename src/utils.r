visium_to_slideseq <- function(sceuratObject){
    ## Description: Change the class of each image to SlideSeq to avoid calling the HE map when drawing
    ## Usage: sceObj <- visium_to_slideseq(sceObj)
    
    for (name in names(sceuratObject@images)){
        sceuratObject@images[[name]] =  new(
                Class = 'SlideSeq',
                assay = "Spatial",
                key = "image_",
                coordinates = sceuratObject@images[[name]]@coordinates[,c('row', 'col')])
    }
    return(sceuratObject)
}
