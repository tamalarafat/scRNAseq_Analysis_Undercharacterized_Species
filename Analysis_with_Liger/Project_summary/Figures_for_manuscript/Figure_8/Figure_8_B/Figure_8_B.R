# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

# Known cell type markers
cell_type_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/known_markers/CH/Curated_markers_hirsuta.csv")

Idents(integrated.data) <- factor(Idents(integrated.data), levels = c(6, 2, 8, 16, 17, 1, 3, 13, 9, 0, 4, 5, 11, 14, 12, 15, 7,10))

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
ggsave(filename = "Figure_8B.png", plot = p, width = 26, height = 20, dpi = 300)

