project_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_markers <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_marker_identification"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_markers, source, .GlobalEnv)

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list_2 <- list.files(paste0(project_dir, "scExplorer/Functions/Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list_2, source, .GlobalEnv)

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

degr = FindMarkers(integrated.data, ident.1 = "7", test.use = "wilcox", only.pos = TRUE)
save(degr, file = "deg_cluster_7.RData")

cdegr = FindConservedMarkers(integrated.data, ident.1 = "7", test.use = "wilcox", grouping.var = "Species", min.pct = 0.1, logfc.threshold = 0.25, only.pos = TRUE)
save(cdegr, file = "cdeg_cluster_7.RData")
