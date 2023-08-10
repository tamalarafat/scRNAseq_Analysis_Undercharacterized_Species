# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library(clevr)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
with = loadRData("withCC.RData")

without = loadRData("withoutCC.RData")

cell_id_with = with$withCC
names(cell_id_with) = with$cell_ID

cell_id_withot = without$withoutCC
names(cell_id_withot) = without$cell_ID

v_measure(cell_id_with, cell_id_withot)
completeness(cell_id_with, cell_id_withot)
homogeneity(cell_id_with, cell_id_withot)
completeness(cell_id_with, cell_id_with)
