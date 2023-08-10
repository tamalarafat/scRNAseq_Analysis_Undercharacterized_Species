# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster 13 - candidate markers
cluster_markers = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_specific_markers/Cluster_13_bio_markers.csv", row.names = 1)

cluster_degs = loadRData("Cluster_13_all_degs.RData")

cluster_degs$significance = "sig"

cluster_degs[!rownames(cluster_degs) %in% rownames(cluster_markers), "significance"] = "nonsig"

cluster_degs$log10 = -log10(cluster_degs$p_val_adj)

cluster_degs[is.infinite(cluster_degs$log10), "log10"] = max(cluster_degs[!is.infinite(cluster_degs$log10), "log10"])

# Load the description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID


# Create ortho genes list and gene description file
description = ATH[rownames(ATH) %in% rownames(cluster_degs), ]

# Description - All DEGs identified for the cluster
markers_description = merge(cluster_degs, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)

markers_description = markers_description[order(markers_description$pct.2, decreasing = FALSE), ]

p <- ggplot(markers_description, aes(x = avg_log2FC, y = log10, color = significance)) +
  geom_point(size = 5) + xlim(-3, 3) + ylim(0, 350) + 
  geom_point(data = markers_description[(markers_description$significance == "sig"), ], aes(x = avg_log2FC, y = log10), size = 10) + 
  scale_color_manual(values = c(grp_col[8], grp_col[3])) + 
  ggrepel::geom_text_repel(
    data = markers_description[(markers_description$significance == "sig") & !is.na(markers_description$TF), ], aes(label = Gene_code_name), nudge_x = 400, size = 16, force = 200, color = "black"
  ) +
  xlab("Average log2FC") + 
  ylab("log10(pval.adj)") + 
  scale_color_manual("Markers", values = c("sig" = grp_col[1], "nonsig" = grp_col[8]), labels = c("Candidates", "Cluster DEGs")) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    axis.title.x = element_text(size = 48, colour = "black"), 
    axis.title.y = element_text(size = 48, colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 48, colour = "black"),
    axis.text.y = element_text(size = 48, colour = "black"),
    title = element_text(size = 40, colour = "black"),
    legend.text = element_text(size = 40, colour = "black"),
    legend.key.size = unit(1.5, 'cm'),
    legend.key.height= unit(0, 'cm'),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position = "top") + 
  guides(color = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "Figure_Volcano_STM.png", plot = p, width = 14, height = 14, dpi = 300)

