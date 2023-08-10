# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library(clevr)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the markers
stm_with_cc = read.csv("Cluster_13_specific_markers_withCC.csv")$gene_ID

stm_without_cc = read.csv("Cluster_13_specific_markers_withoutCC.csv")$gene_ID

length(intersect(stm_with_cc, stm_without_cc))

