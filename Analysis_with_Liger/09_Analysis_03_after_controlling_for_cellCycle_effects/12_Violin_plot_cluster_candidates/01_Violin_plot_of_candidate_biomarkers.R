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


for (i in c(1:length(temp_file_names))){
  
  if (!dir.exists(str_c("Cluster_", parse_number(temp_file_names[i])))){
    dir.create(str_c("Cluster_", parse_number(temp_file_names[i])), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c("Cluster_", parse_number(temp_file_names[i]), "/")
  
  markers = read.csv(str_c(markers_dir, "/", temp_file_names[i]), row.names = 1)
  
  temp_gene_list = split(rownames(markers), ceiling(seq_along(rownames(markers))/9))
  
  for (j in c(1:length(temp_gene_list))){
    
    p <- VlnPlot(integrated.data, features = temp_gene_list[[j]], pt.size = 0)
    
    ggsave(plot = p, filename = str_c(temp_dir, "Cluster_", parse_number(temp_file_names[i]), "_", j, ".png"), width = 12, height = 14, dpi = 300)
    
  }
  
}
