# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library(clevr)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
liger_with = loadRData("withCC.RData")

seurat_with = loadRData("seurat_withCC.RData")

liger_cell_clusters = liger_with$withCC
names(liger_cell_clusters) = liger_with$cell_ID

seurat_cell_clusters = seurat_with$seurat_withCC
names(seurat_cell_clusters) = seurat_with$cell_ID

v_measure(liger_cell_clusters, seurat_cell_clusters)
completeness(liger_cell_clusters, seurat_cell_clusters)
homogeneity(liger_cell_clusters, seurat_cell_clusters)
