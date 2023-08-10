# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

count_table <- table(integrated.data@meta.data$RNA_snn_res.0.3, integrated.data@meta.data$Replicates)
count_mtx   <- as.data.frame.matrix(count_table)
count_mtx$cluster <- rownames(count_mtx)

save(count_mtx, file = "replicate_count_per_cluster.RData")

melt_mtx    <- melt(count_mtx)
melt_mtx$cluster <- as.factor(melt_mtx$cluster)

cluster_size   <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = F))
cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

colnames(melt_mtx)[2] <- "Replicates"

melt_mtx$cluster 

# Set the color for replicates
# base color
base_col = "#f2edee" # The base color

rep_col = colorRampPalette(c(base_col, grp_col[5]))(30)[c(12, 18, 24, 30)]

p <- ggplot(data = melt_mtx, aes(x = cluster, y = value)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = Replicates)) +
  xlab("Clusters") + ylab("Proportion of cells") + 
  coord_flip() + 
  scale_fill_manual(values = c("OX-1" = rep_col[1], "OX-2" = rep_col[2], "OX-3" = rep_col[3], "OX-7" = rep_col[4]), labels = c("Rep 1", "Rep 2", "Rep 3", "Rep 4")) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 48, colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 48, face = "italic")) + 
  guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))

ggsave(filename = "Proportion_replicates.png", plot = p, width = 14, height = 18, dpi = 300)


p <- ggplot(data = melt_mtx, aes(x = cluster, y = value)) + 
  geom_bar(stat = "identity", position = "fill", aes(fill = Replicates)) +
  xlab("Clusters") + ylab("Proportion of cells") + 
  coord_flip() + 
  scale_fill_manual(values = c("OX-1" = grp_col[1], "OX-2" = grp_col[2], "OX-3" = grp_col[5], "OX-7" = grp_col[6]), labels = c("Rep 1", "Rep 2", "Rep 3", "Rep 4")) + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 48, colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"), 
    legend.text = element_text(size = 48, face = "italic")) + 
  guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))

ggsave(filename = "Proportion_replicates_2nd.png", plot = p, width = 14, height = 18, dpi = 300)

