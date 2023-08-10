# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

target_ID = "6"

cluster_cells = WhichCells(integrated.data, idents = target_ID)

# We dont need to load the liger object as the seurat object contains the normalized matrix
inmf_mat = integrated.data@reductions[["inmfcc"]]@cell.embeddings

colnames(inmf_mat) <- str_c("GEP_", parse_number(colnames(inmf_mat)))

inmf_mat = inmf_mat[cluster_cells, ]

# ordered coefficient matrix
df_H <- order_factorized_matrix(inmf_mat)

df_H$GEP <- apply(df_H, 1, function(x) parse_number(colnames(df_H)[which.max(x)]))

table(df_H$GEP)

cluster_GEP_usage = data.frame(table(df_H$GEP))
colnames(cluster_GEP_usage) = c("GEP", "cell_count")

save(cluster_GEP_usage, file = str_c("cluster_", target_ID, "_GEP_usage.RData"))

# Check the STM expressing cells
Cells_with_expression_detection = names(GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", GetAssayData(integrated.data, assay = "RNA", slot = "counts")["Chir06Ox-b35150.2", ] != 0])

writeLines(Cells_with_expression_detection, "RCO_expressing_cells.txt")

cluster_expressing_cells = intersect(Cells_with_expression_detection, cluster_cells)
writeLines(cluster_expressing_cells, "Cluster_6RCO_expressing_cells.txt")

rco_cells_GEP = df_H[cluster_expressing_cells,]
rco_cells_GEP = data.frame(table(rco_cells_GEP$GEP))
save(rco_cells_GEP, file = "rco_cells_GEP.RData")
