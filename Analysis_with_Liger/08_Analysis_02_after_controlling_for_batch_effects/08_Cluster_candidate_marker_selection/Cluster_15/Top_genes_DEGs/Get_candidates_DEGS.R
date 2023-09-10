# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/2023/PhD_projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

# Load the TFs list 
ch_TFs = read.delim("/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_objects/Annotation_files/TFs_list/cardamine_TF_IDs.txt", header = FALSE)[, 1]

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

# Check if the input GEP_genes is a path with the file name. Only .txt file is allowed
DEG_file = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

file_name_pattern = "Cluster"

cluster_IDs = 15


# Create an empty list to store the deg files
temp_deg_list = list()

if (is.character(DEG_file)) {
  
  if (!file_test("-f", DEG_file)) {
    
    # 
    file_names = str_sort(list.files(DEG_file, pattern = file_name_pattern), numeric = TRUE)
    
    for (i in c(1:length(file_names))){
      
      marker_file = loadRData(str_c(DEG_file, "/", file_names[i]))
      
      temp_deg_list[[i]] <- marker_file
      
      if (grepl("[[:digit:]]", file_names[i]) != TRUE){
        names(temp_deg_list)[i] <- sub("(.*)\\.RData$", "\\1", file_names[i])
      }
      
      else {
        names(temp_deg_list)[i] <- str_c("cluster_", parse_number(file_names[i]), "_degs")
      }
    }
  }
  
  else {
    marker_file = loadRData(DEG_file)
    
    temp_deg_list[[1]] <- marker_file
    
    file_names = basename(DEG_file)
    
    if (grepl("[[:digit:]]", file_names) != TRUE){
      names(temp_deg_list)[1] <- sub("(.*)\\.RData$", "\\1", file_names)
    }
    
    else {
      names(temp_deg_list)[i] <- str_c("cluster_", parse_number(file_names), "_degs")
    }
  }
} else {
  if (!is.list(temp_deg_list)){
    temp_deg_list[[1]] <- DEG_file
    
    names(temp_deg_list)[1] <- str_c("cluster_", cluster_IDs, "_degs")
  }
  
  else {
    temp_deg_list = DEG_file
  }
}


cluster_ID = c(15, 16)

find_candidates = 10

# Particular genes set of interest -  could be TFs or other categories
specify_gene_ID = NULL

# Name of the column to be added into the table. If user provides a list of genes, this column will indicate the category -Yes or No 
incorporate_column_name = NULL

# include_both_categories
from_specified_list = TRUE

from_non_specified_list = TRUE

combine_categories = FALSE

# Check if any cluster ID is provided or not. If provided, use the specified ID. If not, get all the cluster IDs and go over each ID.
if (missing(cluster_ID)){
  cluster_ID = levels(Idents(integrated.data))
} else {
  cluster_ID = cluster_ID
}

# Narrow the list of the differentially expressed genes by finding the cluster specific markers
specific_markers_list = specific_marker_finder(DEG_file_dir = DEG_dir, file_name_pattern = "Cluster_", pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE, cluster_ID = cluster_ID)

# If it's a DEG file, get the specific marker for that file

specify_gene_ID = ch_TFs


# Inside the loop now

# What ever the DEG file is put them in a list and go over the list item

# Get the markers file - diffentially expressed genes file
cluster_marker = specific_markers_list[[1]]

# create an empty list to store the files
temp_list = list()

if (!is.null(specify_gene_ID)) {
  
  ## This part for the specified genes set
  
  # First, subset the marker file with the specified gene IDs
  specified_cluster_marker = cluster_marker[cluster_marker$gene_ID %in% specify_gene_ID, ]
  
  # Find candidate markers for the specified genes
  specified_cluster_candidates = select_candidates(specified_cluster_marker)
  
  # Add the source of the markers
  specified_cluster_candidates$source = "DEGs"
  
  # Check if user wants to add the information of the category to the table. If not missing add the column with the given name.
  if (!is.null(incorporate_column_name)) {
    specified_cluster_candidates[[incorporate_column_name]] <- "Yes"
  }
  
  
  ## This part for the remaining genes
  
  # First, subset the marker file with the specified gene IDs
  non_specified_cluster_marker = cluster_marker[!cluster_marker$gene_ID %in% specify_gene_ID, ]
  
  # Find candidate markers for the specified genes
  non_specified_cluster_candidates = select_candidates(non_specified_cluster_marker)
  
  # Add the source of the markers
  specified_cluster_candidates$source = "DEGs"
  
  # Check if user wants to add the information of the category to the table. If not missing add the column with the given name.
  if (!is.null(incorporate_column_name)) {
    non_specified_cluster_candidates[[incorporate_column_name]] <- "NA"
  }
  
  
  ## Now check the requirements of the user.
  
  # If user wants to combine both categories
  if (combine_categories == TRUE) {
    candidate_markers = rbind.data.frame(specified_cluster_candidates, non_specified_cluster_candidates)
    
    # Add the candidates data table to the list
    temp_list[["candidate_markers"]] <- candidate_markers
    
  }
  
  # If user wants to find candidates only from the specified list of genes
  else if (from_specified_list == TRUE) {
    
    # Add the candidates data table to the list
    temp_list[["candidate_markers"]] <- specified_cluster_candidates
    
  }
  
  # If user wants to find candidates only from the specified list of genes
  else if (from_non_specified_list == TRUE) {
    
    # Add the candidates data table to the list
    temp_list[["candidate_markers"]] <- non_specified_cluster_marker
    
  }
} else {
  # Find candidate markers for the specified genes
  candidate_markers = select_candidates(cluster_marker)
  
  # Add the source of the markers
  candidate_markers$source = "DEGs"
  
  # Add the candidates data table to the list
  temp_list[["candidate_markers"]] <- candidate_markers
}

select_candidates <- function(marker_file) {
  
  # Let's define the candidates selection criteria
  marker_file = marker_file[order(marker_file$pct.2, decreasing = FALSE), ]
  
  if (nrow(marker_file) == 0) {
    print("No TFs")} 
  
  else if (nrow(marker_file) <= find_candidates & nrow(marker_file) != 0){
      marker_file = marker_file
    } 
  
  else {
      if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) < find_candidates & nrow(marker_file[marker_file$avg_log2FC >= 1,]) > 0){
        subset1 = marker_file[marker_file$avg_log2FC >= 1,]
        subset2 = marker_file[marker_file$avg_log2FC < 1,]
        subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
        
        if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
          marker_file = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (find_candidates - nrow(subset1)))
          marker_file = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(marker_file[marker_file$avg_log2FC >= 1,]) == 0){
        marker_file = marker_file[order(marker_file$avg_log2FC, decreasing = TRUE), ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
      else {
        marker_file = marker_file[marker_file$avg_log2FC >= 1, ]
        marker_file = marker_file[c(1:find_candidates), ]
      }
    }
  
}




Cluster_marker$source = "DEGs"

write.csv(Cluster_marker, str_c("Cluster_", target_ID, "_top_10_TFs.csv"), row.names = FALSE)


