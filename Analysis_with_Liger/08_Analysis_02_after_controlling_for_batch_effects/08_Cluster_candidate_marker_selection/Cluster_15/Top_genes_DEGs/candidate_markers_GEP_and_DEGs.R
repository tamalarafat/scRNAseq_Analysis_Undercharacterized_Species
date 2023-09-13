projects_dir = "~/Documents/2023/PhD_projects_Yasir/"

# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files(paste0(projects_dir, "scExplorer/Functions"), pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)


# # Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

# load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox")

DEG_files = str_sort(list.files(DEG_dir, pattern = "Cluster"), numeric = TRUE)


# Load the TFs list 
ch_TFs = read.delim(paste0(projects_dir, "Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/cardamine_TF_IDs.txt"), header = FALSE)[, 1]

cluster_IDs = str_sort(levels(integrated.data), numeric = TRUE)

deg_table_list = list()

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


temp_list = list()

for (i in c(1:length(cluster_IDs))) {
  temp_list[[i]] = list(cluster_ID = cluster_IDs[i], 
                        GEP_IDs = GEP_IDs[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(GEP_IDs)))]], 
                        deg_table = deg_table_list[[grep(pattern = str_c("^", cluster_IDs[i], "$", sep = ""), parse_number(names(deg_table_list)))]])
  
  names(temp_list)[i] = str_c("cluster_", cluster_IDs[i])
}

incorporate_column_name = "TFs"

deg_table = loadRData(str_c(DEG_dir, "/", DEG_files[grep(pattern = "^15$", parse_number(DEG_files))]))

candidate_markers_GEP_and_DEGs <- function(seurat_object, 
                                           GEP_IDs = NULL, 
                                           cluster_ID,
                                           DEG_file,
                                           DE_test = "wilcox", 
                                           reduction_name = "inmf",
                                           find_candidates = 10,
                                           max_pct2_detection,
                                           pct_diff,
                                           include_pct_diff,
                                           specify_gene_ID,
                                           incorporate_column_name = NULL,
                                           store_outputs = TRUE,
                                           store_dir = NULL, 
                                           store_folder = "Cluster_candidate_biomarkers") {
  
  
  # Creating necessary storing space to store the results
  
  if (store_outputs == TRUE) {
    # Output directory
    if (missing(store_dir)){
      store_dir = getwd()
    }
    
    # Output folder
    if (!dir.exists(str_c(store_dir, "/", store_folder))){
      dir.create(str_c(store_dir, "/", store_folder), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Output folder for the target cluster
    if (!dir.exists(str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID))){
      dir.create(str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID), showWarnings = TRUE, recursive = FALSE, mode = "0777")
    }
    
    # Assigning the output directory path to a variable
    temp_dir = str_c(store_dir, "/", store_folder, "/", "cluster_", cluster_ID, "/")
  }
  
  if (missing(GEP_IDs)) {
    
    # If GEP ID or IDs are not provided, use the characteristic GEP of the cell cluster to identify candidate markers
    defining_GEP = characteristic_GEP_of_cells(seurat_object = seurat_object, target_ID = cluster_ID, store_output = FALSE)
    
    # Get the ID of the characteristic GEP
    GEP_IDs = as.character(defining_GEP[which.max(defining_GEP$cell_count), "GEP"])
    
  }
  
  # Lets get the 
  gep_candidates = candidate_markers_GEP(seurat_object = seurat_object, 
                                         GEP_IDs = GEP_IDs, 
                                         cluster_ID = cluster_ID, 
                                         find_candidates = find_candidates, 
                                         combine_GEPs = TRUE, 
                                         store_outputs = FALSE)
  
  
  if (!is.null(incorporate_column_name)) {
    gep_candidates[[incorporate_column_name]] = ifelse(gep_candidates$gene_ID %in% specify_gene_ID, "Yes", "No")
  }
  

  
  gep_candidates$source = ifelse(gep_candidates$gene_ID %in% rownames(DEG_file), str_c(gep_candidates$source, "and DEGs", sep = " "), gep_candidates$source)
  
  DEG_file = DEG_file[!rownames(DEG_file) %in% gep_candidates$gene_ID, ]
  
  deg_candidates = candidate_markers_DEGs(DEG_file = DEG_file, 
                                          max_pct2_detection = max_pct2_detection, 
                                          pct_diff = pct_diff, 
                                          include_pct_diff = include_pct_diff, 
                                          find_candidates = find_candidates, 
                                          cluster_ID = cluster_ID, 
                                          specify_gene_ID = specify_gene_ID, 
                                          incorporate_column_name = incorporate_column_name, 
                                          combine_categories = TRUE, 
                                          store_outputs = FALSE)
  
  shared_colnames = colnames(gep_candidates)[colnames(gep_candidates) %in% colnames(deg_candidates)]
  
  gep_candidates = gep_candidates[, shared_colnames]
  
  deg_candidates = deg_candidates[, shared_colnames]
  
  cluster_candidates = rbind.data.frame(gep_candidates, deg_candidates)
  
  # Order the data.frame according to the pct.2 criteria where lowest pct corresponds to the top row / order
  cluster_candidates = cluster_candidates[order(cluster_candidates$pct.2, decreasing = FALSE), ]
  
  # Add rank information - ranks are based on the coefficient. Largest value top rank, lowest value bottom rank
  cluster_candidates$rank = c(1:nrow(cluster_candidates))
  
  return(cluster_candidates)
  
}

cluster_candidates_list = list()

for (i in c(1:length(temp_list))){
  temp_details = temp_list[[i]]
  
  cluster_candidates = candidate_markers_GEP_and_DEGs(seurat_object = integrated.data, 
                                                      GEP_IDs = temp_details$GEP_IDs, 
                                                      cluster_ID = temp_details$cluster_ID, 
                                                      DEG_file = temp_details$deg_table, 
                                                      find_candidates = 10, 
                                                      reduction_name = "inmf", 
                                                      max_pct2_detection = 0.1, 
                                                      pct_diff = 0.3, 
                                                      include_pct_diff = TRUE, 
                                                      specify_gene_ID = ch_TFs, 
                                                      incorporate_column_name = "TFs", 
                                                      store_outputs = FALSE)
  
  cluster_candidates_list[[i]] = cluster_candidates
  
  names(cluster_candidates_list)[i] = str_c(names(temp_list)[i], "candidates", sep = "_")
  
}


cc = candidate_markers_GEP_and_DEGs(seurat_object = integrated.data, GEP_IDs = c(28, 42), cluster_ID = 15, DEG_file = deg_table, find_candidates = 10, reduction_name = "inmf", max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, specify_gene_ID = ch_TFs, incorporate_column_name = "TFs", store_outputs = FALSE)


gep_candidates$source = ifelse(gep_candidates$gene_ID %in% rownames(deg_table), str_c(gep_candidates$source, "and DEGs", sep = " "), gep_candidates$source)

deg_table = deg_table[!rownames(deg_table) %in% gep_candidates$gene_ID, ]

deg_candidates = candidate_markers_DEGs(DEG_file = deg_table, max_pct2_detection = 0.1, pct_diff = 0.3, include_pct_diff = TRUE, find_candidates = 10, cluster_ID = 15, specify_gene_ID = ch_TFs, incorporate_column_name = "TFs", combine_categories = TRUE, store_outputs = FALSE)


shared_colnames = colnames(gep_candidates)[colnames(gep_candidates) %in% colnames(deg_candidates)]

gep_candidates = gep_candidates[, shared_colnames]

deg_candidates = deg_candidates[, shared_colnames]

cluster_candidates = rbind.data.frame(gep_candidates, deg_candidates)



