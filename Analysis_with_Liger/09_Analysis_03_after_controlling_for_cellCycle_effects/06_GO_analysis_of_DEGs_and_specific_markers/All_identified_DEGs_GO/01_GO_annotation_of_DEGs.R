# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

Thaliana_genes =  read.delim("/netscratch/dep_tsiantis/grp_laurent/tamal/2021/Final_Analyses/Input_Files/ATgenes.txt", header = FALSE, col.names = "genes")
Thaliana_genes <- as.character(Thaliana_genes$genes)

storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_GEP20_22_and_CCGEPs"

# load the markers file
load("/home/yasir/Documents/Projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_03_after_controlling_for_cellCycle_effects/Analysis_outputs/DEGs_specific_markers_and_their_orthologues/DEGs_and_Markers.RData")

Cluster_DEGs = markers_list[[2]]

markers_GO_annotation(marker_set = Cluster_DEGs, species_gene_set = Thaliana_genes, store_dir = storing_dir, store_folder = "DEGs_and_Markers_GO", GO_results_folder_name = "All_identified_DEGs_GO")
