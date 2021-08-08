#load data
load('outputData/mtgenes.RData')
load('outputData/mtgenes_distances.RData')

#plot dendrogram
library(cluster)
mtgenes_hca = agnes(dist_ham_genes, diss = T, method = "ward")
dend_mtgenes_hca = as.dendrogram(mtgenes_hca)
par(cex = 0.3) #label size for dendrogram plot
plot(dend_mtgenes_hca, main = "HCA clustering compounds by genes", ylab = "Hamming dist")

#visualize HCA plot
library(ggplot2)
library(factoextra)
fviz_dend(dend_mtgenes_hca,
          k = 5,                     # Cut in 5 groups
          cex = 0.26,                # label size
          k_colors = c("#00AFBB","#E7B800","#9932CC","#458B00","#D02525"),
          color_labels_by_k = TRUE, # color labels by groups 
          main = "HCA clustering compounds by genes",
          ylab = "Hamming dist"
)

#compounds in each group
groups = rect.hclust(mtgenes_hca, k = 5, border = "red")

load('outputData/final_list.RData')
compound_names = names(final_list)

comp_in_group = list()
for (i in 1:5) {
  comp_in_group[[i]] = compound_names[groups[[i]]]
}
save(comp_in_group,file = 'outputData/compounds in HCA groups.RData')













