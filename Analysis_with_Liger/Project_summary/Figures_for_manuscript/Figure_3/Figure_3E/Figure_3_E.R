# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Load the gene file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")


# RCO - AT5G67651 - Chir06Ox-b35150.2
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35150.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_RCO.png", width = 6, height = 6, dpi = 300)


# STM - AT1G62360 - Chir02Ox-b03100.2
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b03100.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_STM.png", width = 6, height = 6, dpi = 300)


# KNAT1 - AT4G08150 - Chir08Ox-b10180.2
p <- FeaturePlot(integrated.data, features = c("Chir08Ox-b10180.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_KNAT1.png", width = 6, height = 6, dpi = 300)


# KNAT2 - AT1G70510 - Chir02Ox-b17040.2
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b17040.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_KNAT2.png", width = 6, height = 6, dpi = 300)


# KNATM - AT1G14760 -Chir01Ox-b14900.2
p <- FeaturePlot(integrated.data, features = c("Chir01Ox-b14900.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_KNATM.png", width = 6, height = 6, dpi = 300)


# LMI1 or ATHB51 - AT5G03790 - Chir06Ox-b35160.2
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35160.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_LMI1.png", width = 6, height = 6, dpi = 300)

# CUC1 - AT3G15170 - Chir03Ox-b14800.2 - ANAC054
p <- FeaturePlot(integrated.data, features = c("Chir03Ox-b14800.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_CUC1.png", width = 6, height = 6, dpi = 300)


# CUC2 - AT5G53950 - Chir08Ox-b15830.2 - ANAC098
p <- FeaturePlot(integrated.data, features = c("Chir08Ox-b15830.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_3E_CUC2.png", width = 6, height = 6, dpi = 300)

