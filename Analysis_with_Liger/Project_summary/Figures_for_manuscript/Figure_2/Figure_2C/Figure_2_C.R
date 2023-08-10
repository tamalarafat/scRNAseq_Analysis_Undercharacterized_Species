# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# get the factor IDs that you would like to plot
factor_ids = c(1:ncol(integrated.data@reductions[["inmf"]]@cell.embeddings))

integrated.data$temp_factor = integrated.data@reductions[["inmf"]]@cell.embeddings[, 20]

integrated.data$temp_factor = (integrated.data$temp_factor)*10

p <- FeaturePlot(integrated.data, features = "temp_factor", reduction = "umap", dims = c(1, 2), pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_2C_GEP20.png", width = 6, height = 6, dpi = 300)


integrated.data$temp_factor = integrated.data@reductions[["inmf"]]@cell.embeddings[, 22]

integrated.data$temp_factor = (integrated.data$temp_factor)*10

p <- FeaturePlot(integrated.data, features = "temp_factor", reduction = "umap", dims = c(1, 2), pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_2C_GEP22.png", width = 6, height = 6, dpi = 300)

