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
res_dir = "/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Liger_analysis/Liger_analysis_with_cell_cycle_effect"

# Create a folder to store/save the output figures
if (!dir.exists(str_c(res_dir, "/", "PCA_plot_replicate_similarity"))){
  dir.create(str_c(res_dir, "/", "PCA_plot_replicate_similarity"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Save the directory information to a variable
temp_dir = str_c(res_dir, "/", "PCA_plot_replicate_similarity", "/")


###
# Load saved data
###
load("/netscratch/dep_tsiantis/common/Yasir/Analysis_of_WT_Cardamine/Liger_analysis/Liger_analysis_with_cell_cycle_effect/Analysis_objects/seurat_object_wt_ox.RData")

# Set the idents to Replicates to generate the PCA plots where cells are colored with the replicate idents
Idents(integrated.data) <- "Replicates"


###
# Type of matrix used to compute PCA - Scaled Raw Data (Default strategy in Seurat)
###

# Log-Normalization of the data
integrated.data <- NormalizeData(integrated.data, normalization.method = "LogNormalize")

# Scale the data
integrated.data <- ScaleData(integrated.data, features = VariableFeatures(integrated.data))

# Run PCA - Compute PCA on the scaled raw data
integrated.data <- RunPCA(integrated.data)

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

ggsave(plot = p, filename = str_c(temp_dir, "PCA_replicate_similarity_scaled_raw_data.png"), width = 14, height = 14, dpi = 300)


###
# Type of matrix used to compute PCA - on the coefficient matrix (reduced dimension of the raw data)
###


# Run PCA - Compute PCA on the coefficient matrix inferred by Liger
integrated.data@reductions["pca"] <- RunPCA(t(integrated.data@reductions$inmfcc@cell.embeddings), npcs = (dim(integrated.data@reductions[["inmfcc"]])[2] - 1))

# Generate plot - coefficient matrix
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

ggsave(plot = p, filename = str_c(temp_dir, "PCA_replicate_similarity_coefficient_matrix.png"), width = 14, height = 14, dpi = 300)


###
# Type of matrix used to compute PCA - on the coefficient matrix (scaled)
###


# To run PCA, prepare the scaled coefficient matrix
integrated.data@reductions[["sinmf"]]@cell.embeddings <- integrated.data@reductions[["sinmf"]]@cell.embeddings[, -c(20, 22)]

# Run PCA - Compute PCA on the scaled coefficient matrix
integrated.data@reductions["pca"] <- RunPCA(t(integrated.data@reductions$sinmf@cell.embeddings), npcs = (dim(integrated.data@reductions[["sinmf"]])[2] - 1))

# Generate plot - scaled coefficient matrix
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

ggsave(plot = p, filename = str_c(temp_dir, "PCA_replicate_similarity_coefficient_matrix_scaled.png"), width = 14, height = 14, dpi = 300)


