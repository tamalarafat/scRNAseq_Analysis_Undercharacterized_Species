# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster GEP
GEP_ID = 23

# get the factor IDs that you would like to plot
factor_ids = colnames(integrated.data@reductions[["inmfcc"]]@cell.embeddings)

factor_index = which(parse_number(colnames(integrated.data@reductions[["inmfcc"]]@cell.embeddings)) == GEP_ID)

integrated.data$temp_factor = integrated.data@reductions[["inmfcc"]]@cell.embeddings[, factor_index]

integrated.data$temp_factor = (integrated.data$temp_factor) * 10

p <- FeaturePlot(integrated.data, features = "temp_factor", reduction = "umap", dims = c(1, 2), pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3I.png", width = 6, height = 6, dpi = 300)
