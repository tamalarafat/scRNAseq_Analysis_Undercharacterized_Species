# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

transgenic_genes = c("Chir01Ox-b14900.2", "Chir02Ox-b16260.2", "Chir03Ox-b00970.2", "Chir03Ox-b16010.2", "Chir03Ox-b17880.2", "Chir04Ox-b19160.2", "Chir05Ox-b03830.2", 
                     "Chir06Ox-b24220.2", "Chir07Ox-b02990.2", "Chir07Ox-b21290.2", "Chir07Ox-b21690.2", "Chir08Ox-b28090.2", "Chir08Ox-b06850.2", "Chir07Ox-b35240.2")

for (i in c(1:length(transgenic_genes))){
  
  p <- FeaturePlot(integrated.data, features = transgenic_genes[i], reduction = "umap", pt.size = 1.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
    NoLegend() + 
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank())
  
  ggsave(plot = p, filename = str_c("Feature_plot_", transgenic_genes[i], ".png"), width = 8, height = 8, dpi = 300)
  
}

# How many genes were picked as cluster markers
# Cluster 5 - candidate markers
cluster_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_specific_markers/Cluster_5_bio_markers.csv", row.names = 1)

shared_genes = intersect(rownames(cluster_markers), transgenic_genes)

for (i in c(1:length(shared_genes))){
  
  p <- FeaturePlot(integrated.data, features = shared_genes[i], reduction = "umap", pt.size = 1.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
    NoLegend() + 
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank())
  
  ggsave(plot = p, filename = str_c("Feature_plot_shared_", shared_genes[i], ".png"), width = 8, height = 8, dpi = 300)
  
}

# Genes
# KNATM = Chir01Ox-b14900.2
# HB13 = Chir02Ox-b16260.2
# HB20 = Chir03Ox-b00970.2
# TIP2-1 = Chir03Ox-b16010.2
# BGLU44 = Chir03Ox-b17880.2
# TIP1-1 = Chir04Ox-b19160.2
# TIP1-2 = Chir05Ox-b03830.2
# HB3 = Chir06Ox-b24220.2
# HB23 = Chir07Ox-b02990.2
# WRKY7 = Chir07Ox-b21290.2
# HB22 = Chir07Ox-b21690.2
# HB25  = Chir08Ox-b28090.2
# CP2 = Chir08Ox-b06850.2
# XTH7 = Chir07Ox-b35240.2
