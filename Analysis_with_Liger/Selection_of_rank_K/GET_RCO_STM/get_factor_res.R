# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# FUnction
get_factor_res <- function(seurat_object_dir, file_name_pattern){
  
  # Lets get the saved file names
  seuratObjects = str_sort(list.files(path = seurat_object_dir, pattern = file_name_pattern), numeric = TRUE)
  
  ml = list()
  
  for (i in c(1:length(seuratObjects))){
    # Creating the file name with path location
    file_name = str_c(seurat_object_dir, "/", seuratObjects[i])
    
    # Loading the liger object and storing it to a temporary variable
    integrated.data = loadRData(file_name) # Load the RData file and assign it to a variable using the function loadRData
    
    DefaultAssay(integrated.data) <- "RNA"
    
    Idents(integrated.data) <- integrated.data$RNA_snn_res.0.3
    
    Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))
    
    integrated.data$cellID = colnames(integrated.data)
    
    # Let's get the cell count plot
    stm_subset = subset(integrated.data, subset = cellID %in% colnames(integrated.data@assays$RNA@counts[, integrated.data@assays$RNA@counts["Chir02Ox-b03100.2", ] != 0]))
    
    stm_cells = names(which.max(table(Idents(stm_subset))))
    
    rco_subset = subset(integrated.data, subset = cellID %in% colnames(integrated.data@assays$RNA@counts[, integrated.data@assays$RNA@counts["Chir06Ox-b35150.2", ] != 0]))
    
    rco_cells = names(which.max(table(Idents(rco_subset))))
    
    
    ml[[i]] <- data.frame(stm_cluster = stm_cells, rco_cluster = rco_cells)
    
    names(ml)[i] <- str_c("K_", dim(integrated.data@reductions[["inmf"]])[2])
  }
  df = do.call(rbind.data.frame, ml)
  
  return(df)
}

# Directory containing the seurat objects
seurat_files_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects"

mdf = get_factor_res(seurat_object_dir = seurat_files_dir, file_name_pattern = "seurat_object_of_")

save(mdf, file = "RCO_STM_Clustering.RData")
