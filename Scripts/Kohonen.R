#Kohonen
load('outputData/mtgenes.RData')
library(kohonen)

#z-value
  ##Should binary data be standardized ???
  #https://www.researchgate.net/post/Should-binary-data-be-standardized/
  #555eb6a36225ffaa518b456b/citation/download. 
  #scale_mt_genes <- scale(mt_genes, center = T, scale = T)

# 1.choosing dimensions by maximum silhouette number ----------------------------------------

source('functions/kohonen_silhouette_function.R')

sil_ta <- grid_size(mt_genes,19,10,inp_distance = "tanimoto")
                                                #max= (1*7)   -0.493372358992498
sil_man <- grid_size(mt_genes,19,10,inp_distance = "manhattan") 
                                                #max= (14*3)  -0.480148887247815
sil_square <- grid_size(mt_genes,19,10,inp_distance = "sumofsquares")
                                                #max= (19*10) -0.63843205328489
sil_eu <- grid_size(mt_genes,19,10,inp_distance = "euclidean") 
                                                #max= (16*9)  -0.653679211887273

# 2.basic supersom --------------------------------------------------------

genes_som <-supersom(mt_genes,grid=somgrid(14,3,"hexagonal"),dist.fcts = "manhattan")


# 3.plotting points -------------------------------------------------------

plot(genes_som, type = "mapping", pchs = 20, main = "genes matrix _ manhattan dist")


# 4.heatmap som -----------------------------------------------------------
colors <- function(n, alpha = 1) rev(heat.colors(n, alpha)) #reverse color ramp 

plot(genes_som, type = "counts", palette.name = colors, heatkey = TRUE,
     main = "genes matrix _ manhattan dist")


#compounds in each cell
load('outputData/final_list.RData')
compound_names = names(final_list)

comp_in_cell = list()
for (i in 1:42) {
  comp_in_cell[[i]] = compound_names[which(genes_som$unit.classif==i)]
}
save(comp_in_cell,file = 'outputData/compounds in kohonen cells.RData')








