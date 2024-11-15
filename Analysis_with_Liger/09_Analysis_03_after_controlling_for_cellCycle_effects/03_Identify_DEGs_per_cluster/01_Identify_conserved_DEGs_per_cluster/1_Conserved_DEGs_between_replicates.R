# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_GEP20_22_and_CCGEPs"

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_GEP20_22_and_CCGEPs/seurat_object_of_K_44_Without_GEP20_and_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

cluster_conserved_marker_finder(seuratObject = integrated.data, grouping_variable = "Genotypes", store_dir = res_dir, DEGtest = "wilcox")
