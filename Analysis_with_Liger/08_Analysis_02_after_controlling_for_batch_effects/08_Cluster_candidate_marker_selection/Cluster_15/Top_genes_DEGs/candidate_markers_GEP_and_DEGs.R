# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/2023/PhD_projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)


# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

cluster_ID = 15

# Lets get the 
gep_candidates = candidate_markers_GEP(seurat_object = integrated.data, GEP_IDs = c(28), cluster_ID = 15, find_candidates = 10, combine_GEPs = TRUE, store_outputs = FALSE)
