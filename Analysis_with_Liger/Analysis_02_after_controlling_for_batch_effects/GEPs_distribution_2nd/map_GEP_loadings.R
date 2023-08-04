# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF")

# Load the seurat object
load("/home/yasir/Documents/Chapter1/Remove_rep_specific_GEP20_and_22/seurat_object_of_K_44_Without_GEP20.RData")


Idents(integrated.data) <- "RNA_snn_res.0.3"

# Lets plot the heatmaps with all the cells
map_factor_loadings(seuratObject = integrated.data, store_dir = getwd(), store_folder = "GEPs_distribution")

