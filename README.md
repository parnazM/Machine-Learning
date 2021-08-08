Unsupervised Machine Learning
===============
The field of machine learning is concerned with the development and application of computer algorithms that improve with experience. Thus, for example, in genomics machine learning can be used to “learn” how to compare and classify compounds which have interactions between genes.  
In Unsupervised approach, the machine learning algorithm takes as input only the unlabeled data and the desired number of different labels to assign.  

Description of the R scripts
---------------
1. preparing binary matrices and compute sil. number   
|R script | binary_matrix.R
|------------ | -------------
|Input data | http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz  ??
|Output data | compounds_genes matrix & compounds_plants matrix  
|            | Silhouette number Table for HCA and Kmeans methods in different cluster numbers and distances
|Package Dependencies | "xlsx" , "e1071"
|Function Dependencies| mesh2name_functon.R , silhouette_functions.R
|Summary | 
---------------
2. visualizing HCA dendrogram   
|R script | HCA.R
|------------ | -------------
|Input data | mtgenes.RData , mtgenes_distances.RData , final_list.RData
|Output data | HCA dendrogram , compounds in each cluster in HCA method
|Package Dependencies | "cluster" , "ggplot2" , "factoextra"
|Summary | 
----------------
3. visualizing Kmeans plot  
|R script | Kmeans.R
|------------ | -------------
|Input data | mtgenes.RData , mtgenes_distances.RData , final_list.RData
|Output data | Kmeans plot , compounds in each cluster in Kmeans method  
|Package Dependencies | "ggplot2" , "factoextra" , "purrr" , "cluster"
|Summary |
