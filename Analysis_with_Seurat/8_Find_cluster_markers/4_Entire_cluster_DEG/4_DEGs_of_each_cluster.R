# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1"

# Load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

cluster_marker_finder(seuratObject = integrated.data, store_dir = res_dir, DEGtest = "wilcox")
