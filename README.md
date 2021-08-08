Unsupervised Machine Learning
===============
The field of machine learning is concerned with the development and application of computer algorithms that improve with experience. Thus, for example, in genomics machine learning can be used to “learn” how to compare and classify compounds which have interactions between genes.  
In Unsupervised approach, the machine learning algorithm takes as input only the unlabeled data and the desired number of different labels to assign.  
Here we compare the genomics similarities of medicinal plants reported to be effective in the management of melancholia depression by using unsupervised machine learning approaches.

Description of the R scripts
---------------
**1. preparing binary matrices and computing sil. number**

|R script | [binary_matrix.R](https://github.com/parnazM/Machine-Learning/blob/master/Scripts/binary_matrix.R)
|------------ | -------------
|Input data | http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz  ??
|Output data | compounds_genes matrix & compounds_plants matrix
|            | Silhouette number Table for HCA and Kmeans methods in different cluster numbers and distances
|Package Dependencies | "xlsx" , "e1071"
|Function Dependencies| mesh2name_functon.R , silhouette_functions.R
|Summary | 
---------------
**2. visualizing HCA dendrogram**   

|R script | [HCA.R](https://github.com/parnazM/Machine-Learning/blob/master/Scripts/HCA.R)
|------------ | -------------
|Input data | mtgenes.RData , mtgenes_distances.RData , final_list.RData
|Output data | HCA dendrogram , compounds in each cluster in HCA method
|Package Dependencies | "cluster" , "ggplot2" , "factoextra"
|Summary | 
----------------
**3. visualizing Kmeans plot**

|R script | [Kmeans.R](https://github.com/parnazM/Machine-Learning/blob/master/Scripts/Kmeans.R)
|------------ | -------------
|Input data | mtgenes.RData , mtgenes_distances.RData , final_list.RData
|Output data | Kmeans plot , compounds in each cluster in Kmeans method  
|Package Dependencies | "ggplot2" , "factoextra" , "purrr" , "cluster"
|Summary |
---------------
**4. finding best distance and dimensions in Kohonen method and visualizing plot**

|R script | [Kohonen.R](https://github.com/parnazM/Machine-Learning/blob/master/Scripts/Kohonen.R)
|------------ | -------------
|Input data | mtgenes.RData , final_list.RData
|Output data | Kohonen plot , heatmap som , compounds in each cluster(cell) in Kohonen method  
|Package Dependencies | "Kohonen" 
|Function Dependencies| kohonen_silhouette_function.R
|Summary |
---------------
**5. trying DBSCAN method for our data**

|R script | [DBSCAN.R](https://github.com/parnazM/Machine-Learning/blob/master/Scripts/DBSCAN.R)
|------------ | -------------
|Input data | mtgenes.RData , mtgenes_distances.RData 
|Output data | DBSCAN plots
|Package Dependencies | "dbscan" , "fpc" , "factoextra"
|Summary |
---------------
