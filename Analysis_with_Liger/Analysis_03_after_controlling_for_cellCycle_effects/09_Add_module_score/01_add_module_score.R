# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Markers directory
markers_dir = "/home/yasir/Documents/Thesis_PhD/Chapter_1/Remove_GEP20_22_and_CCGEPs/Cluster_DEGs_Markers_Candidates/Cluster_biomarkers"

# Lets get the saved file names
temp_file_names = str_sort(list.files(path = markers_dir, pattern = "Cluster"), numeric = TRUE)

marker_list = list()

for (i in c(1:length(temp_file_names))){

  markers = read.csv(str_c(markers_dir, "/", temp_file_names[i]), row.names = 1)
  
  marker_list[[i]] <- rownames(markers)
  
  names(marker_list)[i] <- str_c("cluster", parse_number(temp_file_names[i]), "_module")
}

integrated.data <- AddModuleScore(integrated.data, features = marker_list)

save(integrated.data, file = "/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")


