
# 1.Preparation data such compound & ixns to build binary matrix ----------

# List xls plants
li <- list.files('inputData//thesis/Plants/xls', pattern = '.xls')

# Extract chemical column from plants list
library(xlsx)
h <- list()
for (i in 1:length(li)) {
  f <- read.xlsx2(paste('inputData//thesis//Plants//xls//',li[i],sep = ''), sheetIndex = 1)
  h[[i]] <- f$Chemical
  gc()
}

# Prepare chemical names from plants list
compounds <- gsub(tolower(unique(unlist(h))), pattern = '-', replacement = ' ')

# Extract and prepare chemical names from ixns data
ixns<-read.csv('inputData//CTD_chem_gene_ixns.csv',comment.char = c("#"),
               stringsAsFactors = F) 
colnames(ixns)<-c("ChemicalName","ChemicalID", "CasRN","GeneSymbol", "GeneID","GeneForms",
                  "Organism","OrganismID","Interaction" ,"InteractionActions", "PubMedIDs")
CTD_chemName <- gsub(tolower(unique(ixns$ChemicalName)), pattern = '-', replacement = ' ')

# Intersecte chemical names
intersect_CTD_plants <- intersect(CTD_chemName, compounds)

# Use dictionary 
name_dic <- unique(ixns[,c('ChemicalID','ChemicalName')])
saveRDS(name_dic, file = 'outputData/dictionary.rds')

####
# load('outputData/chem2gene.RData')
#source('functions/mesh2name_functon.R)
# for (i in 1:length(chem2gene)){
#   print(i)
#   names(chem2gene)[i] <- mesh2name(names(chem2gene)[i])
# }
#   
load('outputData/name_chem2gene.RData')
names(chem2gene) <- sapply(names(chem2gene), function(x)gsub(tolower(x), pattern = '-', replacement = ' '))


#Build matrix final_list*number_gene
final_list=chem2gene[intersect_CTD_plants]
final_list<-final_list[-which(final_list=="NULL")]
save(final_list,file='outputData/final_list.RData')

number_genes=unique(unlist(final_list))
mt_genes<-matrix(0, nrow = length(final_list),ncol=length(number_genes))
colnames(mt_genes)=number_genes
rownames(mt_genes)=names(final_list)
for (i in 1:nrow(mt_genes)){
  mt_genes[i,which(colnames(mt_genes) %in% final_list[[i]])]=1
}
  
#Build matrix final_list*plants_name
plants_name<-unlist(li)
plants_name <- gsub(plants_name,pattern = '-farmacy.xls',replacement = '')
plants_name <- gsub(plants_name,pattern = '-',replacement = ' ')
mt_plants<-matrix(0,nrow = length(final_list),ncol = length(plants_name))
rownames(mt_plants)<-names(final_list)
colnames(mt_plants)<-plants_name

for (i in 1:ncol(mt_plants)) {
  sample=h[[i]]
  sample=tolower(sample)
  sample=gsub(sample,pattern = '-',replacement = ' ')
  mt_plants[which(rownames(mt_plants) %in% sample),i]=1
} 
  

save(mt_genes,mt_plants,file='outputData/mtgenes.RData')


# 2. Silhouette number Table for genes matrix -----------------------------

load('outputData/mtgenes.RData')

# Calculate distances in genes matrix
library(cluster)
library(e1071)
dist_eu_genes <- dist(mt_genes, method = "euclidean") #Euclidean distance 
dist_ham_genes <- as.dist(hamming.distance(mt_genes)) #Hamming distance 
dist_ja_genes <- dist(mt_genes, method = "binary")    #Jaccard disance 

save(dist_eu_genes,dist_ham_genes,dist_ja_genes, file = 'outputData/mtgenes_distances.RData')

# Calculate Silouette numb in HCA method for all distance 
source("functions/silhouette_functions.R")
sil_HCA_genes <- lapply(list(dist_eu_genes,dist_ham_genes,dist_ja_genes),
                        function(x)sil_per_clust(x,20,method="HCA"))
sil_HCA_genes <- data.frame(sil_HCA_genes)
colnames(sil_HCA_genes) <- (c('dist_eu_genes','dist_ham_genes','dist_ja_genes'))
save(sil_HCA_genes,file = 'outputData/sil_HCA_genes.RData')

colMax <- function(data) sapply(data, max, na.rm = TRUE) #maximum in each column
colMax(sil_HCA_genes)

# Calculate Silouette numb in KMEANS method for all distance 
sil_KM_genes  <- lapply(list(dist_eu_genes,dist_ham_genes,dist_ja_genes),
              function(x)sil_per_clust(x,20,method="KM"))
sil_KM_genes <- data.frame(sil_KM_genes)
colnames(sil_KM_genes)=(c('dist_eu_genes','dist_ham_genes','dist_ja_genes'))
save(sil_KM,file = 'outputData/sil_KM_genes.RData')

colMax(sil_KM_genes)

table1=cbind(sil_HCA_genes,sil_KM_genes)

library(gplots)
heatmap.2(as.matrix(table1),trace = 'none',srtCol = 45)

# 3. Silhouette number Table for plants matrix ----------------------------

load('outputData/mtgenes.RData')

#calculate distances in plants matrix
dist_eu_plants <- dist(mt_plants, method = "euclidean") #Euclidean distance
dist_ham_plants <- as.dist(hamming.distance(mt_plants)) #Hamming distance
dist_ja_plants <- dist(mt_plants, method = "binary")    #Jaccard distance

save(dist_eu_plants,dist_ja_plants,dist_ham_plants, file = 'outputData/mtplants_distances.RData')

# Calculate Silouette numb in HCA method for all distance
sil_HCA_plants <- lapply(list(dist_eu_plants,dist_ham_plants,dist_ja_plants),
                         function(x)sil_per_clust(x,14,method="HCA"))
sil_HCA_plants <- data.frame(sil_HCA_plants)
colnames(sil_HCA_plants) <- (c('dist_eu_plants','dist_ham_plants','dist_ja_plants'))
save(sil_HCA_plants,file = 'outputData/sil_HCA_plants.RData')

colMax(sil_HCA_plants)

# Calculate Silouette numb in KMEANS method for all distance 
sil_KM_plants=lapply(list(dist_eu_plants,dist_ham_plants,dist_ja_plants),
              function(x)sil_per_clust(x,14,method="KM"))
sil_KM_plants=data.frame(sil_KM_plants)
colnames(sil_KM_plants)=(c('dist_eu_plants','dist_ham_plants','dist_ja_plants'))
save(sil_KM_plants,file = 'outputData/sil_KM_plants.RData')

colMax(sil_KM_plants)

table2=cbind(sil_HCA_plants,sil_KM_plants)

heatmap.2(as.matrix(table2),trace = 'none',srtCol = 45)



# 4.Visualizing matrices------------------------------------

# Clustering genes matrix with ward's method:
##distance by euclidean method
dend_eu_genes <- as.dendrogram(agnes(dist_eu_genes, diss = T, method = "ward"))
plot(dend_eu_genes, main = "clustering compounds by genes", ylab = "Euclidean dist")

##distance by hamming method
dend_ham_genes <- as.dendrogram(agnes(dist_ham_genes, diss = T, method = "ward"))
plot(dend_ham_genes, main = "clustering compounds by genes", ylab = "Hamming dist")

##distance by jaccard method
library(ggplot2)
library(factoextra)
dend_ja_genes <- as.dendrogram(agnes(dist_ja_genes, diss = TRUE, method = "ward"))
fviz_dend(dend_ja_genes,
          k = 3,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c( "#E7B800", "#00AFBB", "#D02525"),
          color_labels_by_k = TRUE, # color labels by groups 
          main = "clustering compounds by genes",
          ylab = "Jaccard dist"
)


# Clustering plants matrix with ward's method
##distance by euclidean method
dend_eu_plants <- as.dendrogram(agnes(dist_eu_plants, diss = T, method = "ward"))
plot(dend_eu_plants, main = "clustering compounds by plants", ylab = "Euclidean dist")

##distance by hamming method
dend_ham_plants <- as.dendrogram(agnes(dist_ham_plants, diss = T, method = "ward"))
plot(dend_ham_plants, main = "clustering compounds by plants", ylab = "Hamming dist")

##distance by jaccard method
dend_ja_plants <- as.dendrogram(agnes(dist_ja_plants, diss = TRUE, method = "ward"))
fviz_dend(dend_ja_plants,
          k = 3,                 # Cut in four groups
          cex = 0.5,                 # label size
          k_colors = c( "#E7B800", "#00AFBB", "#D02525"),
          color_labels_by_k = TRUE, # color labels by groups 
          main = "clustering compounds by plants",
          ylab = "Jaccard dist"
)


# Plot Fan dendrogram for genes matrix
library(ape)
hc_genes <- hclust(dist_eu_genes, method = "ward.D2")
jpeg("plots/dendrograms/genes_fan.jpg", width = 1400, height = 1200, res = 100)
plot(as.phylo(hc_genes), type = "fan", cex = 0.8, main = "Clusterig compounds by genes")
dev.off()

# Plot Fan dendrogram for plants atrix
hc_plants <- hclust(dist_eu_plants, method = "ward.D2")
jpeg("plots/dendrograms/plants_fan.jpg", width = 1400, height = 1200, res = 100)
plot(as.phylo(hc_plants), type = "fan", cex = 0.8, main = "Clusterig compounds by plants")
dev.off()

# Make an Upset plot for plants matrix
library(UpSetR)
df_plants <- as.data.frame(mt_plants)
upset(df_plants, nsets = 14, order.by = "freq", sets.x.label = "Compounds in plants",
      text.scale = c(1.2, 1, 1.2, 1, 1.3, 1.2), shade.color = "red",sets.bar.color = sets_col) 
sets_col <- rev(c("grey51","burlywood3","aquamarine1","dodgerblue2","blue","blueviolet",
                  "gold","darkgoldenrod1","chocolate1","green","deeppink1","maroon",
                  "brown1","brown4")) #color of sets

# 5.Build heatmaps ----------------------------------------------------------------------

library(gplots)
heatmap.2(mt_genes, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
heatmap(mt_genes,scale = "none")

