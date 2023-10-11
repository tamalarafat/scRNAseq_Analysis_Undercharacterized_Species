# projects_dir = "~/Documents/2023/PhD_projects_Yasir/"

projects_dir = "~/Documents/Projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files(paste0(projects_dir, "scExplorer/Functions"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)


# # Load the seurat object
# load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

# Directory path to store the outputs
storing_dir = "/home/yasir/Documents/Projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs"

# Path to the DEG files - differentially genes using the "Findmarkers" function (Seurat) per cluster
DEG_dir = paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox")

# Load the TFs list 
ch_TFs = read.delim(paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/cardamine_TF_IDs.txt"), header = FALSE)[, 1]

# Cluster IDs of the active clustering solution
cluster_IDs = str_sort(levels(integrated.data), numeric = TRUE)

# Empty table to store the degs table
deg_table_list = list()

# From the deg files directory get the file names and load the degs tables
if (is.character(DEG_dir)) {
  
  # get the file names
  file_names = str_sort(list.files(DEG_dir, pattern = "Cluster"), numeric = TRUE)
  
  for (i in c(1:length(file_names))){
    
    marker_file = loadRData(str_c(DEG_dir, "/", file_names[i]))
    
    deg_table_list[[i]] <- marker_file
    
    if (grepl("[[:digit:]]", file_names[i]) != TRUE){
      names(deg_table_list)[i] <- sub("(.*)\\.RData$", "\\1", file_names[i])
    }
    
    else {
      names(deg_table_list)[i] <- str_c("cluster_", parse_number(file_names[i]))
    }
  }
}

# For each cluster ID define the GEP IDs from which the candidate markers will be selected
GEP_IDs = list(cluster_0 = 27, 
               cluster_1 = 5, 
               cluster_2 = 38, 
               cluster_3 = 43, 
               cluster_4 = 44, 
               cluster_5 = 32, 
               cluster_6 = 11, 
               cluster_7 = c(4, 21), 
               cluster_8 = c(29, 33), 
               cluster_9 = 40, 
               cluster_10 = 10, 
               cluster_11 = 34, 
               cluster_12 = 6, 
               cluster_13 = 23, 
               cluster_14 = c(3, 31), 
               cluster_15 = c(28, 42), 
               cluster_16 = 41, 
               cluster_17 = 30)


# Creating an empty list to which GEP IDs and DEGs table for each cluster will be stored. 
# This list will be served as query input to identify candidate biomarkers using the GEP IDs and degs for each cluster
query_list = list()

for (i in c(1:length(cluster_IDs))) {
  query_list[[i]] = list(cluster_ID = cluster_IDs[i], # First element is the cluster ID
                        GEP_IDs = GEP_IDs[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(GEP_IDs)))]], # Second element is the GEP IDs
                        deg_table = deg_table_list[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(deg_table_list)))]]) # 3rd element is the degs table
  
  names(query_list)[i] = str_c("cluster_", cluster_IDs[i])
}


# Now, iterate over the list item to identify candidate biomarkers for each cluster using the query inputs
for (i in c(1:length(query_list))){
  query_inputs = query_list[[i]]
  
  cluster_candidates = candidate_markers_GEP_and_DEGs(seurat_object = integrated.data, 
                                                      GEP_IDs = query_inputs$GEP_IDs, 
                                                      cluster_ID = query_inputs$cluster_ID, 
                                                      DEG_file = query_inputs$deg_table, 
                                                      find_candidates = 10, 
                                                      reduction_name = "inmf", 
                                                      max_pct2_detection = 0.1, 
                                                      pct_diff = 0.3, 
                                                      include_pct_diff = TRUE, 
                                                      specify_gene_ID = ch_TFs, 
                                                      incorporate_column_name = "TFs", 
                                                      store_outputs = TRUE, 
                                                      store_dir = storing_dir)
}
