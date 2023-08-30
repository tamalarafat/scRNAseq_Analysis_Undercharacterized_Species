# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/Yasirs_projects/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("~/Documents/Beginning_of_a_compendium/Chapter_2/comparative_study_of_two_species/Liger_Analysis_strategy_2/Seurat_object/seurat_object_of_K_50.RData")

# Get the characteristic GEP of cluster 5
cluster_GEP_usage = characteristic_GEP_of_cells(seurat_object = integrated.data, 
                                                gene_ID =  "Chir06Ox-b35150.2", 
                                                reduction_name = "inmf", 
                                                store_output = TRUE, 
                                                store_folder = "Defining_GEPs_of_RCO_cells")

# Max GEP of RCO TF expressing cells
paste0("Defining GEP of RCO TF expressing cells: GEP ", as.character(cluster_GEP_usage[which.max(cluster_GEP_usage$cell_count), "GEP"]))
