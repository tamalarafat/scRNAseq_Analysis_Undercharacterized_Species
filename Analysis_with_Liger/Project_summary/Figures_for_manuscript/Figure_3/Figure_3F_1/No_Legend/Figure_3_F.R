# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster 13 cells
cluster_cells = WhichCells(integrated.data, idents = "13")

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_cells, ]

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

axis_text = as.character(c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43))

axis_text_face = rep("plain", length(axis_text))

axis_text_face[which(axis_text == "23")] = "bold"

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
  scale_y_discrete(breaks = str_c("GEP_", c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_3_STM_cluster.png", plot = p, width = 12, height = 18, dpi = 300)

######
# 2nd heatmap without sorting
######
# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cluster_cells, ]

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
  scale_y_discrete(breaks = str_c("GEP_", c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_3_STM_cluster_without_sorting.png", plot = p, width = 12, height = 18, dpi = 300)


######
# 3rd heatmap - rco cells sorted
######

# Check the STM expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir02Ox-b03100.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir02Ox-b03100.2", ] != 0])

writeLines(Cells_with_expression_detection, "STM_expressing_cells.txt")

cluster_expressing_cells = intersect(Cells_with_expression_detection, cluster_cells)

writeLines(cluster_expressing_cells, "STM_expressing_cells_cluster_13.txt")

non_expressing_cells = setdiff(cluster_cells, Cells_with_expression_detection)

cells_ordered = c(cluster_expressing_cells, non_expressing_cells)

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
  scale_y_discrete(breaks = str_c("GEP_", c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black"))  + NoLegend()

ggsave(filename = "Figure_3_STM_cluster_stm_cells_sorted.png", plot = p, width = 12, height = 18, dpi = 300) + NoLegend()


######
# 4th heatmap - only rco expressing cells
######
# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H = df_H[cluster_expressing_cells, ]

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
  xlab("Cluster 13 - stm expressing cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43)), labels = axis_text) +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 42, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 42, colour = "black", face = axis_text_face),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 42, colour = "black")) + NoLegend()

ggsave(filename = "Figure_3_STM_cluster_only_stm_cells.png", plot = p, width = 12, height = 18, dpi = 300) + NoLegend()

