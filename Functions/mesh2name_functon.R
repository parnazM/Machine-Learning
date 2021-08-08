
# Function that take chemical ID , give chemical Name

mesh2name<-function(inp){
   dictionary <- readRDS('outputData/dictionary.rds')
   final_name=''
   for (i in 1:nrow(dictionary)){
   if (inp==dictionary$ChemicalID[i]){
    final_name=dictionary$ChemicalName[i]}  
   } 
  return(final_name)
}