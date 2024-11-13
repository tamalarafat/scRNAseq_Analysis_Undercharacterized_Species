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


# Output directory
res_dir = "/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Seurat_analysis/Replicate_similarity_plots"

# Create a folder to store/save the output figures
if (!dir.exists(str_c(res_dir, "/", "PCA_plot_replicate_similarity"))){
  dir.create(str_c(res_dir, "/", "PCA_plot_replicate_similarity"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Save the directory information to a variable
temp_dir = str_c(res_dir, "/", "PCA_plot_replicate_similarity", "/")


###
# Load saved data
###
load("/netscratch/dep_tsiantis/common/scRNAseq/yasir/WT_Cardamine/Seurat_Analysis/Integrated_seurat_object/integrated_ox_wt_seurat.RData")

integrated.data$Replicates <- stri_sub(integrated.data$Replicates, 4, -1)

integrated.data$Replicates <- as.factor(integrated.data$Replicates)

Idents(integrated.data) <- "Replicates"

###
# Default strategy in Seurat
###


# Generate plot - scaled raw data
p <- DimPlot(integrated.data, reduction = "pca", dims = c(1, 2), label = F, pt.size = 1.5, repel = T, label.size = 20) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"),
    legend.key.size = unit(4,"line"),
    legend.key.height = unit(0.5, 'cm'),
    legend.key.width = unit(0.5, 'cm'),
    legend.key = element_rect(size = 20),
    legend.spacing = unit(0.05, 'cm'),
    legend.text = element_text(size = 18, face = "bold"),
    legend.spacing.x = unit(0.1, 'cm')) +
  guides(colour = guide_legend(override.aes = list(size = 8)))

ggsave(plot = p, filename = str_c(temp_dir, "PCA_replicate_similarity.png"), width = 14, height = 14, dpi = 300)
