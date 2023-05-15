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


rotate_2D_coordinate <- function(coordinate_mat, angle){
    ## Description: Rotate 2d coordinates　https://math.stackexchange.com/questions/195141/rotating-2d-coordinates
    ## param coordinate_mat: x,y coordinate matrix
    ## param angle: angle
    ## Usage: mat <- rotate_2D_coordinate(coordinate_mat = mat, angle = 30)
    
    # Calculate the sine and cosine of the angle
    sin_a <- sin(angle * pi/180)
    cos_a <- cos(angle * pi/180)
    # Define the rotation matrix
    rot_mat <- matrix(c(cos_a, -sin_a, sin_a, cos_a), nrow=2, ncol=2, byrow=TRUE)   
    # Rotate the array using the rotation matrix
    rot_a <- t(sapply(1:nrow(coordinate_mat), function(i){rot_mat  %*% coordinate_mat[i,]})) 
    return (rot_a)
}


## 10x Visium images object
# 10X HE相对于表达矩阵是颠倒的，如果做伪图，需要注意。另外imagerow:Y, imagecol:X
# 表达矩阵x-y保持不变，改变image的矩阵，具体
# 1.imagerow:Y, imagecol:X
# 2.对imagerow做颠倒
# coord$x_r <- max(coord$x) - coord$x + min(coord$x)


if (FALSE){
DefaultImage <- function(object) {
  object <- UpdateSlots(object = object)
  images <- Images(object = object, assay = DefaultAssay(object = object))
  if (length(x = images) < 1) {
    images <- Images(object = object)
  }
  return(images[[1]])
}

    
image <- DefaultImage(seu)
seu@images[[image]]@coordinates <- seu@meta.data[,c('y', 'x')]
colnames(seu@images[[image]]@coordinates) <- c('imagerow', 'imagecol')
seu@images[[image]]@coordinates$imagerow <- max(seu@images[[image]]@coordinates$imagerow) - seu@images[[image]]@coordinates$imagerow + min(seu@images[['slice1']]@coordinates$imagerow)

    
    
p1 <- SpatialFeaturePlot(seu, features =  'Spp1', pt.size.factor = 1, 
                   images = NULL,  stroke = NA, alpha = c(1, 1)) + 
    theme_cowplot()+
    coord_fixed()


p2 <- ggplot(seu@meta.data, aes(x = x, y =y))+
    geom_point(shape = 19, size = 0.1)+
    theme_cowplot()+
    coord_fixed()

}
