project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Storing directory
res_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1"

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

# Convert the class of the grouping variable - character to factor
integrated.data$Genotype <- factor(integrated.data$Genotype)

# Set the cluster idents of the cells
Idents(integrated.data) <- "integrated_snn_res.0.4"

# Identify conserved markers for the replicates
cluster_conserved_marker_finder(seuratObject = integrated.data, grouping_variable = "Genotype", store_dir = res_dir, DEGtest = "wilcox", between_conditions = FALSE)
