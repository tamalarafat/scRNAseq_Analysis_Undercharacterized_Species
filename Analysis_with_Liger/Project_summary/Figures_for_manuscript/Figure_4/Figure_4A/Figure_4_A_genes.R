# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster 5 - candidate markers
cluster_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_specific_markers/Cluster_5_bio_markers.csv", row.names = 1)

# Create a markers list to plot
feature_markers = c("Chir06Ox-b35150.2", "Chir08Ox-b07590.2", "Chir06Ox-b35870.2", "Chir05Ox-b16230.2", "Chir01Ox-b35290.2", "Chir03Ox-b30470.2", 
       "Chir04Ox-b20280.2", "Chir01Ox-b29650.2", "Chir08Ox-b06850.2", "Chir01Ox-b42060.2", "Chir07Ox-b11520.2")

markers = cluster_markers[feature_markers, ]

for (j in c(1:nrow(markers))){
    p <- FeaturePlot(integrated.data, features = rownames(markers)[j], reduction = "umap", pt.size = 1.5, order = T, min.cutoff = 0.001, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])) + 
      NoLegend() + 
      theme(
        line = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        plot.title = element_blank())
    
    ggsave(plot = p, filename = if (!is.na(markers[j, "Gene_code_name"])){str_c("Feature_plot_", markers[j, "Gene_code_name"], ".png")} else {str_c("Feature_plot_", rownames(markers)[j], ".png")}, width = 6, height = 6, dpi = 300)
    }

