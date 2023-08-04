# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/home/yasir/Documents/Chapter1/Remove_rep_specific_GEP20_and_22/seurat_object_of_K_44_Without_GEP20.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# read the markers file
markers = read.csv("Cluster_7_bio_markers.csv", row.names = 1)

for (i in c(1:nrow(markers))){
  
  if (!dir.exists("markers_expression")){
    dir.create("Markers_expression", showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  p <- FeaturePlot(integrated.data, features = rownames(markers)[i], reduction = "umap", pt.size = 1.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
    NoLegend() + 
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank())
  
  ggsave(plot = p, filename = if (!is.na(markers[i, "Gene_code_name"])){str_c("Markers_expression/", "Feature_plot_", markers[i, "Gene_code_name"], ".png")} else {str_c("Markers_expression/", "Feature_plot_", rownames(markers)[i], ".png")}, width = 8, height = 8, dpi = 300)
  
}

