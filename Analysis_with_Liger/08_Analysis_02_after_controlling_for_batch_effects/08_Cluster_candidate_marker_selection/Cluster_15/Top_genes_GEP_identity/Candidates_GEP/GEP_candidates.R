# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/Yasirs_projects/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("~/Documents/Beginning_of_a_compendium/Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Liger_objects_with_CC/seurat_object_of_K_44_Without_GEP20_22.RData")

# Lets get the 
candidate_markers_GEP(seurat_object = integrated.data, GEP_IDs = c(28, 42), cluster_ID = 15, find_candidates = 10, combine_GEPs = TRUE)

