# Load all the functions stored in scripts from the folder housing the scripts
library(ggplot2)
library(Seurat)
## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"


# RCO - AT5G67651 - Chir06Ox-b35150.2
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35150.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_RCO.png", width = 6, height = 6, dpi = 300)


# STM - AT1G62360 - Chir02Ox-b03100.2
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b03100.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_STM.png", width = 6, height = 6, dpi = 300)


# KNAT1 - AT4G08150 - Chir08Ox-b10180.2
p <- FeaturePlot(integrated.data, features = c("Chir08Ox-b10180.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_KNAT1.png", width = 6, height = 6, dpi = 300)


# KNAT2 - AT1G70510 - Chir02Ox-b17040.2
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b17040.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_KNAT2.png", width = 6, height = 6, dpi = 300)


# KNATM - AT1G14760 -Chir01Ox-b14900.2
p <- FeaturePlot(integrated.data, features = c("Chir01Ox-b14900.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_KNATM.png", width = 6, height = 6, dpi = 300)


# LMI1 or ATHB51 - AT5G03790 - Chir06Ox-b35160.2
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35160.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_LMI1.png", width = 6, height = 6, dpi = 300)

# CUC1 - AT3G15170 - Chir03Ox-b14800.2 - ANAC054
p <- FeaturePlot(integrated.data, features = c("Chir03Ox-b14800.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_CUC1.png", width = 6, height = 6, dpi = 300)


# CUC2 - AT5G53950 - Chir08Ox-b15830.2 - ANAC098
p <- FeaturePlot(integrated.data, features = c("Chir08Ox-b15830.2"), reduction = "umap", pt.size = 2, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
  NoLegend() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    plot.title = element_blank())

ggsave(plot = p, filename = "Figure_8F_CUC2.png", width = 6, height = 6, dpi = 300)

