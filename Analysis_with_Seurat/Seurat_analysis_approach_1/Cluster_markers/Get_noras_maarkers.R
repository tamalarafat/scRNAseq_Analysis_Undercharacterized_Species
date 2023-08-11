# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Dir - containing the DEG files
DEG_file = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/Markers_Nora/Chirsuta_WT_leaf_lognorm_cca_yasir0.4_wilcox_markers.csv", row.names = 1)

# RCO
Cluster_DEGs = DEG_file[DEG_file$cluster == "7", ]

# Prepare the file of cluster differentially expressed gene
Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
rownames(Cluster_DEGs) = Cluster_DEGs$Chirsutav12
Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)

# Here create a function for specific marker selection
Cluster_specific_DEGs = Cluster_DEGs[(Cluster_DEGs$pct.2 <= 0.1) | ((Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= 0.3), ]
Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)

# Save the file
write.csv(Cluster_specific_DEGs, file = "seurat_rco_cluster_7.csv", row.names = FALSE)


# STM
Cluster_DEGs = DEG_file[DEG_file$cluster == "15", ]

# Prepare the file of cluster differentially expressed gene
Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
rownames(Cluster_DEGs) = Cluster_DEGs$Chirsutav12
Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)

# Here create a function for specific marker selection
Cluster_specific_DEGs = Cluster_DEGs[(Cluster_DEGs$pct.2 <= 0.1) | ((Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= 0.3), ]
Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)

# Save the file
write.csv(Cluster_specific_DEGs, file = "seurat_stm_cluster_15.csv", row.names = FALSE)
