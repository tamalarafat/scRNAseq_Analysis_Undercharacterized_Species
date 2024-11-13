# Specify the projects dir containing all the functions of scRNA-seq analysis pipeline
projects_dir <- "/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions/"

# Functions - Markers identification
list_1 <- list.files(paste0(projects_dir, "Library_handler"), pattern = "*.R$", full.names = TRUE)
sapply(list_1, source, .GlobalEnv)

# Functions - Data manipulation
list_2 <- list.files(paste0(projects_dir, "Functions_matrix_manipulation"), pattern = "*.R$", full.names = TRUE)
sapply(list_2, source, .GlobalEnv)

# Functions - Library and packages handler
list_3 <- list.files(paste0(projects_dir, "Functions_marker_identification"), pattern = "*.R$", full.names = TRUE)
sapply(list_3, source, .GlobalEnv)

## Define the colors
grp_col = c("#00A08A", "#075149FF")

# Output directory
res_dir = "/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Seurat_analysis/Replicate_similarity_plots"

# Create a folder to store/save the output figures
if (!dir.exists(str_c(res_dir, "/", "Proportion_plot_replicate_similarity"))){
  dir.create(str_c(res_dir, "/", "Proportion_plot_replicate_similarity"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Save the directory information to a variable
temp_dir = str_c(res_dir, "/", "Proportion_plot_replicate_similarity", "/")


###
# Load saved data
###
load("/netscratch/dep_tsiantis/common/scRNAseq/yasir/WT_Cardamine/Seurat_Analysis/Integrated_seurat_object/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

integrated.data$Replicates <- stri_sub(integrated.data$Replicates, 4, -1)

integrated.data$Replicates <- as.factor(integrated.data$Replicates)

count_table <- table(integrated.data@meta.data$integrated_snn_res.0.4, integrated.data@meta.data$Replicates)
count_mtx   <- as.data.frame.matrix(count_table)
count_mtx$cluster <- rownames(count_mtx)

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

rep_col = colorRampPalette(c(base_col, grp_col[2]))(30)[c(12, 18, 24, 30)]

p <- ggplot(data = melt_mtx, aes(x = cluster, y = value)) +
  geom_bar(stat = "identity", position = "fill", aes(fill = Replicates)) +
  xlab("Clusters") + ylab("Proportion of cells") +
  coord_flip() +
  scale_fill_manual(values = c("OX-1" = rep_col[1], "OX-2" = rep_col[2], "OX-3" = rep_col[3], "OX-7" = rep_col[4]), labels = c("Rep 1", "Rep 2", "Rep 3", "Rep 4")) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[1], size = 2),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.text = element_text(size = 48, colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 48, face = "italic")) +
  guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))

ggsave(filename = str_c(temp_dir, "Proportion_replicates.png"), plot = p, width = 14, height = 18, dpi = 300)


p <- ggplot(data = melt_mtx, aes(x = cluster, y = value)) +
  geom_bar(stat = "identity", position = "fill", aes(fill = Replicates)) +
  xlab("Clusters") + ylab("Proportion of cells") +
  coord_flip() +
  scale_fill_manual(values = c("OX-1" = rep_col[1], "OX-2" = rep_col[2], "OX-3" = rep_col[3], "OX-7" = rep_col[4])) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = grp_col[1], size = 2),
    panel.background = element_blank(), # Background of the entire plot
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.text = element_text(size = 48, colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(2, "line"),
    legend.text = element_text(size = 48, face = "italic")) +
  guides(fill = guide_legend(title = "Species"), color = guide_legend(override.aes = list(size = 8)))

ggsave(filename = str_c(temp_dir, "Proportion_replicates_exp.png"), plot = p, width = 14, height = 18, dpi = 300)

