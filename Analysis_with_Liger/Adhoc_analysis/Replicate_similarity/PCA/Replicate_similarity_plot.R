# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/"

# Functions - Markers identification
list_1 <- list.files(paste0(projects_dir, "Library_handler"), pattern = "*.R$", full.names = TRUE) 
sapply(list_1, source, .GlobalEnv)

# Functions - Data manipulation
list_2 <- list.files(paste0(projects_dir, "Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(list_2, source, .GlobalEnv)

# Functions - Library and packages handler
list_3 <- list.files(paste0(projects_dir, "Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(list_3, source, .GlobalEnv)

###
# Load saved data
###

load("/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Liger_analysis/Liger_analysis_with_cell_cycle_effect/Analysis_objects/seurat_object_wt_ox.RData")

# Set the idents to Replicates to generate the PCA plots where cells are colored with the replicate idents
Idents(integrated.data) <- "Replicates"

# Output directory
res_dir = "/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Liger_analysis/Liger_analysis_with_cell_cycle_effect"

# Generate the PCA plot
RDimension_plot(seuratObject = integrated.data, store_dir = res_dir, store_folder = "PCA_plot_replicate_similarity", dimension_reduction_name = "pca")
