project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/Conserved_DEGs_grouped_by_Replicates/Conserved_markers_DEtest_wilcox")

# Storing directory
storing_dir = paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs")

# Generate and store the specific markers for the cell clusters
specific_marker_finder(DEG_file_dir = DEG_dir, 
                       file_name_pattern = "Cluster", 
                       max_pct2_detection = 0.1, 
                       pct_diff = 0.3, 
                       include_pct_diff = TRUE, 
                       store_folder = "Cluster_specific_conserved_markers", 
                       store_dir = storing_dir, 
                       store_outputs = TRUE)
