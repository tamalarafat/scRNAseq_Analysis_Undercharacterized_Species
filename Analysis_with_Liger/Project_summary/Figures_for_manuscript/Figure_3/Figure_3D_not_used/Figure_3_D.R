# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

#### Lets have a visualization of species view
Idents(integrated.data) <- "Replicates"

# base color
base_col = "#f2edee" # The base color

rep_col = colorRampPalette(c(base_col, grp_col[5]))(30)[c(12, 18, 24, 30)]

# No legend
p <- DimPlot(integrated.data, reduction = "umap", dims = c(1, 2), label = F, pt.size = 2, repel = T, shuffle = TRUE) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 42, face = "italic")) + 
  scale_color_manual(values = c("OX-1" = grp_col[1], "OX-2" = grp_col[2], "OX-3" = grp_col[5], "OX-7" = grp_col[6]), labels = c("Rep 1", "Rep 2", "Rep 3", "Rep 4")) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Figure_3D.png", width = 12, height = 12, dpi = 300)


