# Get GEPs from the generated Basis file - containing grouping of the genes
candidate_markers_DEGs <- function(seurat_object, 
                                   GEP_genes = NULL, 
                                   GEP_IDs = NULL, 
                                   cluster_ID,
                                   DE_test = "wilcox", 
                                   reduction_name = "inmf",
                                   find_candidates = 10,
                                   GEP_source_name = NULL,
                                   combine_GEPs = FALSE,
                                   store_outputs = TRUE,
                                   store_dir = NULL, 
                                   store_folder = "DEG_candidates_of_clusters"){
  
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
  
  
  
  # Get the genes of the GEP of interest. User can specify the genes as
  # 1 - Without specifying any GEP IDs or set of genes. If not specified, the characteristic GEP of the cluster will be picked to get the genes
  # 2 - File path including the text file. Only .txt is allowed
  # 3 - A vector containing the gene IDs
  # 4 - GEP ID or IDs
  
  # This part will take the argument input and get the genes of the GEP or GEPs and assign them into a list that will be used later to iterate over the GEP IDs.
  
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
  
  
  # create an empty list to store the files
  temp_list = list()
  
  
  # Check if user specified a set of genes to look for candidates from the DEGs file
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
  } 
  
  else {
    # Find candidate markers for the specified genes
    candidate_markers = select_candidates(cluster_marker)
    
    # Add the source of the markers
    candidate_markers$source = "DEGs"
    
    # Add the candidates data table to the list
    temp_list[["candidate_markers"]] <- candidate_markers
  }
  

}

if (store_outputs == TRUE) {
  
}
