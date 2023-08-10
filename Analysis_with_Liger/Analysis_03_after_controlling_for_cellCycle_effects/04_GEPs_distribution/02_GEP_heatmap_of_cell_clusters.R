# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"
Idents(integrated.data) <- factor(Idents(integrated.data), levels = seq(0, length(levels(Idents(integrated.data))) - 1))

cluster.ids = str_sort(levels(Idents(integrated.data)), numeric = TRUE)

for (i in c(1:length(cluster.ids))){
  heatmap_cells_coefficient(seuratObject = integrated.data, store_dir = getwd(), reduction_pattern = "inmfcc",store_folder = "Heatmap_of_clusters_inmfcc", cell_ids = WhichCells(integrated.data, idents = cluster.ids[i]), figureName = str_c("Cluster_", cluster.ids[i]))
}

save(integrated.data, file = "/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

md = integrated.data@meta.data
md$cell_ID = rownames(md)

md = md[, c("cell_ID", "RNA_snn_res.0.3")]
colnames(md)[2] <- "withoutCC"

save(md, file = "withoutCC.RData")
