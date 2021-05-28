
#Silhouette bar plot function for 1:k
sil_per_clust <- function(input_distance, k, method="HCA"){
  require(purrr)
  require(cluster)
  
  sil=rep(NA,k)  
  
  if(method == "HCA") {
    sil[1]=0
    for (i in 2:k) {
      trees <- cutree(agnes(input_distance, diss = T, method = "ward"), k = i)
      sil[i]=mean(silhouette(trees, dist = input_distance)[,3])
    } #end for
    
  } else {
    sil[1]=0
    for (i in 2:k) {
      km <- kmeans(input_distance, centers = i, nstart = 50)
      sil[i]=mean(silhouette(km$cluster, input_distance)[,3])
      
    }
  } #end for
  
  plot(sil, type = "l", xlab = "K", main = paste(method,"method"))
  abline(v=i[which.max(sil)], col="red",  lty=1)    #????????????
  return(sil)
  
} # end function


# best iter.max number for 1:k , kmeans method
km_iter.max=function(input_distance=dist_ham_genes, k=10 ,max_iter=1){
  require(cluster)
  
  sil_width=matrix(NA,nrow = k,ncol = max_iter)
  
  for (i in 2:k){
    for (j in 1:max_iter){
      model_km <- kmeans(input_distance, centers = i,iter.max = j+10, nstart = 50)
      sil_width[i,j]=mean(silhouette(model_km$cluster, input_distance)[,3])
      
      
    } # end for j
  }#end for i
  
  names(sil_width[i,]) = paste(c("k",k[i]))   # ????????
  names(sil_width[,j]) = max_iter[j+10]      # ??????????
  
  w <- which(sil_width == max(sil_width, na.rm = T), arr.ind = TRUE)
  print(paste("k=", w[,1], "and iter.max= ", (w[,2]+10), "have most silhouette number !"))
  return(sil_width)
  
}# end function
