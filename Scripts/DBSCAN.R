#DBSCAN
load('outputData/mtgenes.RData') #load matrix data
load('outputData/mtgenes_distances.RData')

  ##Large value of EPS , may form broad clusters , reduce the number of clusters
   #and may result into reduction in number of noise point .

  ##Small value of EPS may introduce more noise points.

  ##Large value of MinNumPoints introduce more noise but will form robust clusters.

  ##Small value of MinNumPoints includes noise poitns in clusters , border points
   #can be misclassified & considered as core points..

# 1. optimum value of minpoints -------------------------------------------

  #For two-dimensional data: use default value of minPts=4 (Ester et al., 1996)

# 2. optimum value of eps -------------------------------------------------

library(dbscan)
kNNdistplot(dist_ham_genes, k=4)      #k = Minpoints
abline(h=70, lty=2, col= "red")      #first knee => eps = 70

kNNdistplot(dist_ja_genes, k=4)      #k = Minpoints
abline(h=0.005, lty=2, col= "red")  #eps = 0.005

kNNdistplot(dist_eu_genes, k=4)      #k = Minpoints
abline(h=0.1, lty=2, col= "red")    #eps = 0.1

kNNdistplot(mt_plants, k=4)
abline(h=0.01, lty=2, col= "red") # first knee => eps = 0.01


# 3. Vizualize DBSCAN plots ------------------------------------------------------------

# compute DBSCAN 
library(fpc)
db_ham_genes <- dbscan(dist_ham_genes,
                  eps= 70,         #Reachability maximum distance
                  MinPts = 4,        #Reachability minimum number of points
                  scale = FALSE, 
                  method = "dist")

db_ja_genes <- dbscan(dist_ja_genes,
                   eps= 0.005,         #Reachability maximum distance
                   MinPts = 4,        #Reachability minimum number of points
                   scale = FALSE, 
                   method = "dist")

db_eu_genes <- dbscan(dist_eu_genes,
                   eps= 0.1,         #Reachability maximum distance
                   MinPts = 4,        #Reachability minimum number of points
                   scale = FALSE, 
                   method = "dist")


db_plants <- dbscan(mt_plants,
                   eps= 0.05,        #Reachability maximum distance
                   MinPts = 4,       #Reachability minimum number of points
                   scale = FALSE, 
                   method = "dist")

# Plot DBSCAN results using factoextra package
library(factoextra)
fviz_cluster(db_eu_genes, mt_genes, stand = FALSE, ellipse = TRUE, geom = "point",
             main = 'DBSCAN _ euclidean dist')
fviz_cluster(db_ham_genes, mt_genes, stand = FALSE, ellipse = TRUE, geom = "point",
             main = 'DBSCAN _ hamming dist')
fviz_cluster(db_ja_genes, mt_genes, stand = FALSE, ellipse = TRUE, geom = "point",
             main = 'DBSCAN _ jaccard dist')

