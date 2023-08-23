features_weight_across_clusters <- function(features_list, store_dir = NULL, store_folder = "Duplicated_genes_details") {
  
  # Creating necessary storing space to store the results
  
  if (missing(store_dir)){
    store_dir = getwd()
  }
  
  if (!dir.exists(str_c(store_dir, "/", store_folder))){
    dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
  }
  
  temp_dir = str_c(store_dir, "/", store_folder, "/")
  
  if(features_list == "character") {
    
  }
}


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

