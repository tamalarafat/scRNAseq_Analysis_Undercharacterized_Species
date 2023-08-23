# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

library(readxl)
library(english)
library(Seurat)


# Load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

# read the xlsx file

genes_file = as.data.frame(read_excel(path = "List_duplicates_Cardamine.xlsx", col_names = "Chr_genes"))

genes_file$gene_ID = str_replace(string = sub(pattern = "^(.*_.*)_.*", replacement = "\\1", x = genes_file$Chr_genes), pattern = "_", replacement = "-")

df_present = genes_file[genes_file$gene_ID %in% intersect(genes_file$gene_ID, rownames(integrated.data)), ]

df_present[, c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = ""

rownames(df_present) = df_present$gene_ID

df_absent = genes_file[!genes_file$gene_ID %in% intersect(genes_file$gene_ID, rownames(integrated.data)), ]

df_absent[, c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = ""

# Lets create a empty vector to store the count information
temp = list()

cluster_labels = str_sort(levels(Idents(integrated.data)), numeric = TRUE)

search_max = 2

for (i in c(1:length(cluster_labels))){
  cell_names = WhichCells(integrated.data, idents = cluster_labels[i])
  cell_count = sum(Matrix(GetAssayData(integrated.data, assay = "RNA", slot = "counts"), sparse = TRUE)[df_present$gene_ID[1], cell_names] != 0)
  temp[[i]] = data.frame(cluster_ID = (i-1), expressing_cells = cell_count)
}

df_count_holder = do.call(rbind.data.frame, temp)

df_count_holder = df_count_holder[order(df_count_holder$expressing_cells, decreasing = TRUE)[1:search_max], ]

df_present[df_present$gene_ID[1], c(str_c("cluster_ID_max_", english(c(1:search_max))), str_c("cell_count_max_", english(c(1:search_max))))] = c(df_count_holder$cluster_ID[c(1:search_max)], df_count_holder$expressing_cells[c(1:search_max)])

