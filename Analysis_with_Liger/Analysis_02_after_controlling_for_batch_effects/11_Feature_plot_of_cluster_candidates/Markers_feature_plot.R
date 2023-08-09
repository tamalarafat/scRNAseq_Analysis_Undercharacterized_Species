# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/home/yasir/Documents/Chapter1/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Markers directory
markers_dir = "/home/yasir/Documents/Chapter1/Liger_analysis_strategy_1/Results/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_specific_markers"

# Lets get the saved file names
temp_file_names = str_sort(list.files(path = markers_dir, pattern = "Cluster"), numeric = TRUE)

# Loading the liger object and storing it to a temporary variable

for (i in c(1:length(temp_file_names))){
  if (!dir.exists(str_c("Cluster_", parse_number(temp_file_names[i])))){
    dir.create(str_c("Cluster_", parse_number(temp_file_names[i])), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c("Cluster_", parse_number(temp_file_names[i]), "/")
  
  markers = read.csv(str_c(markers_dir, "/", temp_file_names[i]), row.names = 1)
  
  for (j in c(1:nrow(markers))){
    p <- FeaturePlot(integrated.data, features = rownames(markers)[j], reduction = "umap", pt.size = 1.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
      NoLegend() + 
      theme(
        line = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank())
    
    ggsave(plot = p, filename = if (!is.na(markers[j, "Gene_code_name"])){str_c(temp_dir, "Feature_plot_", markers[j, "Gene_code_name"], ".png")} else {str_c(temp_dir, "Feature_plot_", rownames(markers)[j], ".png")}, width = 8, height = 8, dpi = 300)
    
  }
  
}
