# kohonen dimensions function ----------------------------------------------------------------

grid_size <- function(inp_matrix,dim1,dim2,inp_distance = "tanimoto"){
  require(kohonen)
  require(cluster)
  
  sil <- matrix(NA,nrow = dim1, ncol = dim2)
  for (i in 1:dim1) {
    for (j in 2:dim2) {
      som_dist <-supersom(inp_matrix, grid=somgrid(xdim = i,ydim = j,"hexagonal"),
                          dist.fcts = inp_distance)
      sil_width <- silhouette(som_dist$unit.classif, dist = som_dist$distances)[,3]
      sil_width <- sil_width[!(is.infinite(sil_width) | is.nan(sil_width))]
      #print(sil_width)
      #print(length(which(is.nan(sil_width) | is.infinite(sil_width))))
      sil[i,j] <- mean(sil_width,na.rm = T)
      print(paste(i,j,sil[i,j]))
      
    }#end for
  }#end for
  
  return(sil)
  
} #end function


