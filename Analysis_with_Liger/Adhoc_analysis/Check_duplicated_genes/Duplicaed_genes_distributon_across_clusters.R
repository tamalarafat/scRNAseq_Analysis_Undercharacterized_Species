library(readxl)
library(english)
library(Seurat)
library(stringr)


# Load the seurat object
load("Your seurat object")

# read the xlsx file

input_file = as.data.frame(read_excel(path = "List_duplicates_Cardamine.xlsx", col_names = "Chr_genes"))

input_file$gene_ID = str_replace(string = sub(pattern = "^(.*_.*)_.*", replacement = "\\1", x = input_file$Chr_genes), pattern = "_", replacement = "-")

# Call the function with desired parameters
features_weight_across_clusters(seuratObject = integrated.data, features_list = a, search_max = 2, include_cell_count = TRUE, store_dir = getwd())
