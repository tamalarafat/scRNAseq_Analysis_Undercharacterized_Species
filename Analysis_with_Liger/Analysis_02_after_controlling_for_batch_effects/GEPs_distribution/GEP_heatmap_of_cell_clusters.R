# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

cluster.ids = levels(Idents(integrated.data))


for (i in c(1:length(cluster.ids))){
  heatmap_cells_coefficient(seuratObject = integrated.data, store_dir = getwd(), reduction_pattern = "inmfcc",store_folder = "Heatmap_of_clusters_inmfcc", cell_ids = WhichCells(integrated.data, idents = cluster.ids[i]), figureName = str_c("Cluster_", cluster.ids[i]))
}

md = integrated.data@meta.data
md$cell_ID = rownames(md)

md = md[, c("cell_ID", "RNA_snn_res.0.3")]
colnames(md)[2] <- "withCC"

save(md, file = "withCC.RData")
