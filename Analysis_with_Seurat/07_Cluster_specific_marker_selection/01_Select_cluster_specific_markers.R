project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Seurat/Analysis_outputs/Differentially_expressed_genes/DEG_DEtest_wilcox")

# Storing directory
storing_dir = paste0(project_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Seurat/Analysis_outputs")

specific_marker_finder(DEG_file_dir = DEG_dir, max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, store_dir = storing_dir, store_outputs = TRUE)

# Load the seurat object
sob = readRDS("/netscratch/dep_tsiantis/common/scRNAseq/160123_Newest_AthCol0-ChWT-rco_res_from_Nora/generated.Rdata/integrated_ox_wt_seurat.rds")

load("/netscratch/dep_tsiantis/common/scRNAseq/yasir/WT_Cardamine/Seurat_Analysis/Integrated_seurat_object/integrated_ox_wt_seurat.RData")
