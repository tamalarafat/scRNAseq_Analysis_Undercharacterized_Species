# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# With legend
p <- DimPlot(integrated.data, reduction = "umap", dims = c(1, 2), label = T, pt.size = 2, repel = T, label.size = 20) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Figure_2A.png", width = 14, height = 14, dpi = 300)


# No legend
p <- DimPlot(integrated.data, reduction = "umap", dims = c(1, 2), label = T, pt.size = 2, repel = T, label.size = 20) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank()
    )

ggsave(plot = p, filename = "Figure_2A_2nd.png", width = 14, height = 14, dpi = 300)
