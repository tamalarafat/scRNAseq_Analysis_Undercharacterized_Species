# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Known cell type markers
cell_type_markers = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/known_markers/CH/Curated_markers_hirsuta.csv")

Idents(integrated.data) <- factor(Idents(integrated.data), levels = c(0, 8, 15, 16, 1, 2, 3, 4, 6, 9, 10, 17, 12, 11, 14, 13, 5, 7))

# Dotplot
p <- DotPlot(integrated.data, features = cell_type_markers$CH_ID, scale = T, assay = "RNA", dot.scale = 25, col.min = 0, dot.min = 0.05, cols = RColorBrewer::brewer.pal(9, name = "Greens")[c(2,9)]) +
  guides(size = guide_legend(override.aes = list(color= RColorBrewer::brewer.pal(9, name = "Greens")[9], alpha = 1))) +
  xlab("Genes") +
  ylab("Cell clusters") +
  scale_x_discrete(labels = as.character(cell_type_markers$Gene_Name)) + 
  theme_classic() + 
  theme(axis.title.x = element_text(size = 42),
        axis.title.y = element_text(size = 42),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(color = grp_col[2], size = 2),
        axis.text.x = element_text(size = 42, colour = "black", angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 42, colour = "black"),
        legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 34, colour = "black"),
        legend.spacing = unit(3.0, 'cm'), 
        legend.position = "top", 
        legend.key.width = unit(3.0, 'cm'),
        legend.title = element_text(size = 42, colour = "black"),
        legend.box.spacing = unit(2.5, 'cm'))
ggsave(filename = "Figure_3C.png", plot = p, width = 26, height = 20, dpi = 300)

