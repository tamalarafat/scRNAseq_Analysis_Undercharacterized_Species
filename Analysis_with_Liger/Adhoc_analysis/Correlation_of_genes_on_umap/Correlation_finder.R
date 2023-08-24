library(Seurat)

# Load the Seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

# For gene coexpression analysis, please refer to this site:: https://github.com/satijalab/seurat/issues/519

# CUC3 / NAC / ANAC031 - Chir02Ox-b23540.2 - AT1G76420
# For CUC3, find coexpression to Chir01Ox-b42590.2 (AT1G53700 - WAG1) and Chir03Ox-b14070.2 (AT3G14370 - WAG2)

# CUC3 - WAG1
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b23540.2", "Chir01Ox-b42590.2"), blend = TRUE, reduction = "umap", order = TRUE, pt.size = 1) + 
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

ggsave(plot = p, filename = "Feature_plot_CUC3_WAG1.png", width = 28, height = 12, dpi = 300)


# CUC3 - WAG2
p <- FeaturePlot(integrated.data, features = c("Chir02Ox-b23540.2", "Chir03Ox-b14070.2"), blend = TRUE, reduction = "umap", order = TRUE, pt.size = 1) + 
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

ggsave(plot = p, filename = "Feature_plot_CUC3_WAG2.png", width = 28, height = 12, dpi = 300)



# RCO - Chir06Ox-b35150.2 - AT5G67651
# For RCO, find correlation to Chir04Ox-b06770.2 (AT2G24260 - LRL1), Chir06Ox-b08920.2 (AT5G48170 - SNE), Chir07Ox-b10470.2 (AT4G14550 - SLR)

# RCO - LRL1
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35150.2", "Chir04Ox-b06770.2"), blend = TRUE, reduction = "umap", order = TRUE, pt.size = 1) + 
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

ggsave(plot = p, filename = "Feature_plot_RCO_LRL1.png", width = 28, height = 12, dpi = 300)


# RCO - SNE
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35150.2", "Chir06Ox-b08920.2"), blend = TRUE, reduction = "umap", order = TRUE, pt.size = 1) + 
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

ggsave(plot = p, filename = "Feature_plot_RCO_SNE.png", width = 28, height = 12, dpi = 300)


# RCO - SLR
p <- FeaturePlot(integrated.data, features = c("Chir06Ox-b35150.2", "Chir07Ox-b10470.2"), blend = TRUE, reduction = "umap", order = TRUE, pt.size = 1) + 
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

ggsave(plot = p, filename = "Feature_plot_RCO_SLR.png", width = 28, height = 12, dpi = 300)

####
# Feature scatter for correlation coefficient
####

#### This part show the correlation values of the genes

Idents(integrated.data) <- integrated.data$Replicates

p <- FeatureScatter(integrated.data, feature1 = "Chir06Ox-b35150.2", feature2 = "Chir04Ox-b06770.2", slot = "data", cols = NULL) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Correlation_plot_RCO_LRL1.png", width = 12, height = 12, dpi = 300, bg = "white")

p <- FeatureScatter(integrated.data, feature1 = "Chir06Ox-b35150.2", feature2 = SNE, slot = "data", cols = NULL) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Correlation_plot_RCO_SNE.png", width = 12, height = 12, dpi = 300, bg = "white")

p <- FeatureScatter(integrated.data, feature1 = "Chir06Ox-b35150.2", feature2 = SLR, slot = "data", cols = NULL) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Correlation_plot_RCO_SLR.png", width = 12, height = 12, dpi = 300, bg = "white")

p <- FeatureScatter(integrated.data, feature1 = CUC3, feature2 = WAG1, slot = "data", cols = NULL) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Correlation_plot_CUC3_WAG1.png", width = 12, height = 12, dpi = 300, bg = "white")

p <- FeatureScatter(integrated.data, feature1 = CUC3, feature2 = WAG2, slot = "data", cols = NULL) + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) + 
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = "Correlation_plot_CUC3_WAG2.png", width = 12, height = 12, dpi = 300, bg = "white")


# Genes
RCO = "Chir06Ox-b35150.2"
LRL1 = "Chir04Ox-b06770.2"
SNE = "Chir06Ox-b08920.2"
SLR = "Chir07Ox-b10470.2"
CUC3 = "Chir02Ox-b23540.2" 
WAG1 = "Chir01Ox-b42590.2"
WAG2 = "Chir03Ox-b14070.2"

# Get the data from the Seurat object
matrix <- as.matrix(GetAssayData(integrated.data, slot = "counts", assay = "RNA"))

cor(matrix[RCO, ], matrix[LRL1, ]) # -0.009767335
cor(matrix[RCO, ], matrix[SNE, ]) # -0.001492406
cor(matrix[RCO, ], matrix[SLR, ]) # -0.006468677
cor(matrix[CUC3, ], matrix[WAG1, ]) # 0.1277713
cor(matrix[CUC3, ], matrix[WAG2, ]) # -0.003855503


