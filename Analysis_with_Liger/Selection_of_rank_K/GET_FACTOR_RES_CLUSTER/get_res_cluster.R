# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# FUnction
get_res_cluster <- function(seurat_object_dir, file_name_pattern){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  # This list to store seurat object results
  temp_list = list()
  
  for (j in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[j])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    DefaultAssay(integrated.data) <- "RNA"
    
    # let's go over the resolution values
    temp = str_sort(names(integrated.data@meta.data)[str_detect(names(integrated.data@meta.data), pattern = "res")], numeric = TRUE)
    
    # this list is to save resolution results for each seurat objet
    ml = list()
    
    for (i in c(1:length(temp))){
      
      temp_str = substr(gsub("[^[:alnum:]]", "", temp[i]), start = nchar(gsub("[^[:alnum:]]", "", temp[i])) - 1, stop = nchar(gsub("[^[:alnum:]]", "", temp[i])))
      temp_str = str_c("RES", gsub(pattern = "[A-Za-z]", replacement = "", temp_str))
      
      # Set the cluster identity to the seurat object
      Idents(integrated.data) <- temp[i]
      
      ml[[i]] <- length(levels(Idents(integrated.data)))
      
      names(ml)[i] <- temp_str
    }
    
    df = do.call(cbind.data.frame, ml)
    
    temp_list[[j]] <- df
    
    names(temp_list)[j] <- str_c("K_", dim(integrated.data@reductions[["inmf"]])[2])
  }
  
  res_df = do.call(rbind.data.frame, temp_list)
  
  return(res_df)
}

# Directory containing the seurat objects
seurat_files_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects"

mdf = get_res_cluster(seurat_object_dir = seurat_files_dir, file_name_pattern = "seurat_object_of_")

save(mdf, file = "Factor_res_clusters.RData")


