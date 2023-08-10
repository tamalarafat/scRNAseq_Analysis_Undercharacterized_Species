# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library("ggVennDiagram")

library(VennDiagram)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

liger_rco = read.csv("Cluster_5_specific_markers_description.csv")

liger_stm = read.csv("Cluster_13_specific_markers_description.csv")

seurat_rco = read.csv("seurat_rco_cluster_7.csv")

seurat_stm = read.csv("seurat_stm_cluster_15.csv")


x = list(liger_rco = liger_rco$gene_ID,
                   liger_stm = liger_stm$gene_ID,
                   seurat_rco = seurat_rco$Chirsutav12,
                   seurat_stm = seurat_stm$Chirsutav12)

p <- ggvenn(
  x, 
  fill_color = c(grp_col[2], grp_col[1], grp_col[4], grp_col[9]),
  fill_alpha = 0.8,
  stroke_size = 0.8, 
  set_name_size = 0,
  text_size = 8
)

ggsave(plot = p, filename = "Figure_8E.png", width = 10, height = 10, dpi = 300, bg = "white")

p <- ggvenn(
  x, 
  fill_color = c(grp_col[2], grp_col[1], grp_col[4], grp_col[9]),
  fill_alpha = 0.8,
  stroke_size = 0.8, 
  set_name_size = 5,
  text_size = 8
)

ggsave(plot = p, filename = "Figure_8E_2.png", width = 10, height = 10, dpi = 300, bg = "white")
