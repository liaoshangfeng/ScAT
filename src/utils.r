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

## rotate 2d coordinates　https://math.stackexchange.com/questions/195141/rotating-2d-coordinates
rotate_2D_coordinate <- function(coordinate_mat, angle){
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
