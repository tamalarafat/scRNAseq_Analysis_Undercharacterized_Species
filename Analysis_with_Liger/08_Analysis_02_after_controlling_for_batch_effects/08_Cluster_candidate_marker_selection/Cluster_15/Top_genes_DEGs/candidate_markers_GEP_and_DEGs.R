projects_dir = "~/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files(paste0(projects_dir, "scExplorer/Functions"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)


# # Load the seurat object
# load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox")

# Load the TFs list 
ch_TFs = read.delim(paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/cardamine_TF_IDs.txt"), header = FALSE)[, 1]

cluster_ID = 15

if (!is.null(incorporate_column_name)) {
  gep_candidates[[incorporate_column_name]] = "Yes/No" - based on the condition
}
  
# Lets get the 
gep_candidates = candidate_markers_GEP(seurat_object = integrated.data, GEP_IDs = c(28), cluster_ID = 15, find_candidates = 10, combine_GEPs = TRUE, store_outputs = FALSE)

deg_candidates = candidate_markers_DEGs(DEG_file = DEG_dir, file_name_pattern = "Cluster", max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, find_candidates = 10, cluster_ID = 15, specify_gene_ID = ch_TFs, incorporate_column_name = "TFs", combine_categories = TRUE, store_outputs = FALSE)

shared_colnames = colnames(gep_candidates)[colnames(gep_candidates) %in% colnames(deg_candidates)]