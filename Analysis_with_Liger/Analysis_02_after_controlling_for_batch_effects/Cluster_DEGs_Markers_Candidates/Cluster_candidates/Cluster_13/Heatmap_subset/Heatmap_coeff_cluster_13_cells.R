# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library(CountClust)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", hue_pal()(47)[20], "#767676FF", "#FD8CC1FF")

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster 8 cells
cluster_13_cells = WhichCells(integrated.data, idents = "13")

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_13_cells, ]

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% dplyr::select(-GEP) %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 13 cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(23, 44)), labels = c(23, 44)) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black"),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, face = "bold"))

ggsave(filename = "Heatmap_with_cluster_13_GEPs.png", plot = p, width = 12, height = 18, dpi = 300)


######
# 2nd heatmap without sorting
######
# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cluster_13_cells, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 13 cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(23, 44)), labels = c(23, 44)) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black"),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, face = "bold"))

ggsave(filename = "Heatmap_with_cluster_13_GEPs_not_sorted.png", plot = p, width = 12, height = 18, dpi = 300)


# Check the STM expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir02Ox-b03100.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir02Ox-b03100.2", ] != 0])

writeLines(Cells_with_expression_detection, "STM_expressing_cells.txt")

STM_expressing_cells = intersect(Cells_with_expression_detection, cluster_13_cells)

rest_cells <- setdiff(cluster_13_cells, STM_expressing_cells)

cells_ordered = c(STM_expressing_cells, rest_cells)

######
# 3rd heatmap without sorting
######
# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[STM_expressing_cells, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 13 cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(23, 44)), labels = c(23, 44)) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black"),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, face = "bold"))

ggsave(filename = "Heatmap_with_cluster_13_GEPs_stm_cells_sorted.png", plot = p, width = 12, height = 18, dpi = 300)


######
# 3rd heatmap without sorting
######
# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cells_ordered, ]

df_H$cell_id = rownames(df_H)
df_H$cell_id = factor(df_H$cell_id, levels = df_H$cell_id)

df_H <-
  df_H %>% melt(
    id.vars = "cell_id",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(df_H, aes(x = cell_id, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion of\nGEP usage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("Cluster 13 cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(23, 44)), labels = c(23, 44)) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, face = "bold", color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black"),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, face = "bold", colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, face = "bold"))

ggsave(filename = "Heatmap_with_cluster_13_GEPs_only_stm_cells.png", plot = p, width = 12, height = 18, dpi = 300)
