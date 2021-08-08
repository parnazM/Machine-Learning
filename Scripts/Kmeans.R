# K_Mean clustering
load('outputData/mtgenes.RData')
load('outputData/mtgenes_distances.RData')

# visualize kmeans cluster ------------------------------------------------

km_model = kmeans(dist_ham_genes,centers = 2, nstart = 50)

library(ggplot2)
library(factoextra)
fviz_cluster(km_model, data = dist_ham_genes,
             palette = c("#00AFBB", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

#compounds in each cluster
comp_in_clust = list(
names(which(km_model$cluster==1)),
names(which(km_model$cluster==2))
)
save(comp_in_clust,file = 'outputData/compounds in kmeans clusters.RData')

# elbow plot --------------------------------------------------------------
library(purrr)
library(ggplot2)

tot_withinss<- map_dbl(1:10, function(k){
  model <- kmeans(x = dist_ham_genes, centers = k)
  model$tot.withinss
})
elbow_df_genes <- data.frame(k=1:10, tot_withinss=tot_withinss)
ggplot(elbow_df_genes, aes(x = k, y = tot_withinss)) +
  geom_line() +
  scale_x_continuous(breaks = 1:10)


# Silhouette analysis -----------------------------------------------------

library(purrr)
library(cluster)
library(ggplot2)

sil_width2 <- map_dbl(2:10, function(k){
  model <- pam(mt_genes, k = k)
  model$silinfo$avg.width
})
sil_df <- data.frame(
  k = 2:10,
  sil_width2 = sil_width2
)
print(sil_df)
ggplot(sil_df, aes(x=k, y=sil_width2))+
  geom_line()+
  scale_x_continuous(breaks = 2:10)
