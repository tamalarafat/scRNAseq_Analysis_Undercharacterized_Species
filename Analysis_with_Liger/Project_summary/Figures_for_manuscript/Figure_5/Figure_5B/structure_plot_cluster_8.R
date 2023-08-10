# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

library(CountClust)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", hue_pal()(47)[20], "#767676FF", "#FD8CC1FF")

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Cluster 8 cells
cluster_8_cells = WhichCells(integrated.data, idents = "8")

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_8_cells, ]


#######
# Structure plot 1 -  All GEPs
#######

mdf <- order_factorized_matrix(inmf_mat)
mdf = cbind(mdf[, str_c("GEP_", c(5, 11, 27, 29, 33, 44))], "GEP_sum" = rowSums(mdf[, !colnames(mdf) %in% str_c("GEP_", c(5, 11, 27, 29, 33, 44))]))
mdf$GEP <- apply(mdf, 1, function(x) colnames(mdf)[which.max(x)])


cell_clusters = mdf[, "GEP", drop = FALSE]
cell_clusters$sample_id = rownames(cell_clusters)
cell_clusters = cell_clusters[, c(2, 1)]
cell_clusters <- cell_clusters[str_order(cell_clusters$GEP, decreasing = T, numeric = TRUE),]
colnames(cell_clusters)[2] = "tissue_label"
rownames(cell_clusters) <- NULL

cell_clusters$tissue_label <- as.factor(cell_clusters$tissue_label)

mat = mdf[cell_clusters$sample_id, -8]

# Structure plot
p1 <- StructureGGplot(
  omega = mat,
  palette = grp_col[c(9, 3, 1, 2, 5, 6, 8)],
  annotation = cell_clusters,
  yaxis_label = "Cell Clusters",
  order_sample = TRUE,
  sample_order_decreasing = T,
  legend_title_size = 16,
  legend_key_size = 0.6,
  legend_text_size = 16
) +
  ylab("Membership proportion") +
  xlab("Cell clusters") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,face = "bold", colour = "black"), 
        axis.ticks.length = unit(0, "cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())
ggsave(plot = p1, filename = str_c("SP_CL8_GEP_", paste(c(5, 11, 27, 29, 33, 44), collapse = "_"), ".png"), width = 10, height = 12)



#######
# Structure plot 2 - without 5
#######

mdf <- order_factorized_matrix(inmf_mat)
mdf = cbind(mdf[, str_c("GEP_", c(11, 27, 29, 33, 44))], "GEP_sum" = rowSums(mdf[, !colnames(mdf) %in% str_c("GEP_", c(11, 27, 29, 33, 44))]))
mdf$GEP <- apply(mdf, 1, function(x) colnames(mdf)[which.max(x)])


cell_clusters = mdf[, "GEP", drop = FALSE]
cell_clusters$sample_id = rownames(cell_clusters)
cell_clusters = cell_clusters[, c(2, 1)]
cell_clusters <- cell_clusters[str_order(cell_clusters$GEP, decreasing = T, numeric = TRUE),]
colnames(cell_clusters)[2] = "tissue_label"
rownames(cell_clusters) <- NULL

cell_clusters$tissue_label <- as.factor(cell_clusters$tissue_label)

mat = mdf[cell_clusters$sample_id, -7]

# Structure plot
p1 <- StructureGGplot(
  omega = mat,
  palette = grp_col[c(3, 1, 2, 5, 6, 8)],
  annotation = cell_clusters,
  yaxis_label = "Cell Clusters",
  order_sample = TRUE,
  sample_order_decreasing = T,
  legend_title_size = 16,
  legend_key_size = 0.6,
  legend_text_size = 16
) +
  ylab("Membership proportion") +
  xlab("Cell clusters") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 16,face = "bold", colour = "black"), 
        axis.ticks.length = unit(0, "cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.background = element_blank())
ggsave(plot = p1, filename = str_c("SP_CL8_GEP_", paste(c(11, 27, 29, 33, 44), collapse = "_"), ".png"), width = 10, height = 12)



#######
# Structure plot 3 - active GEPs
#######
mdf <- order_factorized_matrix(inmf_mat)
mdf = mdf[, str_c("GEP_", c(11, 27, 29, 33, 44))]
mdf <- order_factorized_matrix(mdf)
mdf$GEP <- apply(mdf, 1, function(x) colnames(mdf)[which.max(x)])


cell_clusters = mdf[, "GEP", drop = FALSE]
cell_clusters$sample_id = rownames(cell_clusters)
cell_clusters = cell_clusters[, c(2, 1)]
cell_clusters <- cell_clusters[str_order(cell_clusters$GEP, decreasing = T, numeric = TRUE),]
colnames(cell_clusters)[2] = "tissue_label"
rownames(cell_clusters) <- NULL

cell_clusters$tissue_label <- as.factor(cell_clusters$tissue_label)

mat = mdf[cell_clusters$sample_id, -6]

# Structure plot
p1 <- StructureGGplot(
  omega = mat,
  palette = grp_col[c(3, 1, 2, 5, 6)],
  annotation = cell_clusters,
  yaxis_label = "Cell Clusters",
  order_sample = TRUE,
  sample_order_decreasing = T,
  legend_title_size = 16,
  legend_key_size = 0.6,
  legend_text_size = 16) +
  scale_fill_manual(name = "GEPs", values = c("GEP_11" = grp_col[3], "GEP_27" = grp_col[1], "GEP_29" = grp_col[2], "GEP_33" = grp_col[5], "GEP_44" = grp_col[6]), 
                      labels = str_c("GEP ", c(11, 27, 29, 33, 44))) + 
  ylab("Membership proportion") +
  xlab("Cell clusters") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32, colour = "black"), 
        axis.ticks.length = unit(0, "cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 32, colour = "black"),
        legend.title = element_text(size = 32, colour = "black"),
        panel.background = element_blank()) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 4)))
ggsave(plot = p1, filename = str_c("SP_CL8_GEP_only_", paste(c(11, 27, 29, 33, 44), collapse = "_"), ".png"), width = 10, height = 12)


# No legend
# Structure plot
p1 <- StructureGGplot(
  omega = mat,
  palette = grp_col[c(3, 1, 2, 5, 6)],
  annotation = cell_clusters,
  yaxis_label = "Cell Clusters",
  order_sample = TRUE,
  sample_order_decreasing = T,
  legend_title_size = 16,
  legend_key_size = 0.6,
  legend_text_size = 16
) + NoLegend() + 
  labs(fill = "GEPs") + 
  ylab("Membership proportion") +
  xlab("Cell clusters") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 32, colour = "black"), 
        axis.ticks.length = unit(0, "cm"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 32, colour = "black"),
        legend.title = element_text(size = 32, colour = "black"),
        panel.background = element_blank())
ggsave(plot = p1, filename = str_c("No_legend_SP_CL8_GEP_only_", paste(c(11, 27, 29, 33, 44), collapse = "_"), ".png"), width = 10, height = 12)

