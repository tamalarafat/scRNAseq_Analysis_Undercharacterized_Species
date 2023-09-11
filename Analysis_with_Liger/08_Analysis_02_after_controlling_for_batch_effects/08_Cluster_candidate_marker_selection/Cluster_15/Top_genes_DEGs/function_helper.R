# Get GEPs from the generated Basis file - containing grouping of the genes
candidate_markers_DEGs <- function(DEG_file,
                                   file_name_pattern = "Cluster",
                                   pct_detection = 0.1,
                                   pct_diff = 0.3,
                                   find_candidates = 10,
                                   cluster_ID = NULL,
                                   specify_gene_ID = NULL,
                                   incorporate_column_name = NULL,
                                   from_specified_list = FALSE,
                                   from_non_specified_list = FALSE,
                                   combine_categories = FALSE,
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
  
  # A nested function to select candidates from the markers file
  select_candidates <- function(marker_file, find_candidates = 10) {
    
    n_pct = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "pct.2")]
    
    if (!"gene_ID" %in% colnames(cluster_degs)){
      marker_file$gene_ID = rownames(marker_file)
    }
    
    
    if (length(n_pct) == 1){
      # Let's define the candidates selection criteria
      marker_file = marker_file[order(marker_file[["pct.2"]], decreasing = FALSE), ]
      
      if (nrow(marker_file) == 0) {
        print("No TFs")
      } 
      
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
    
    else if (length(n_pct) == 2){
      marker_file = marker_file[order(marker_file[, n_pct[1]], marker_file[, n_pct[2]]), ]
      rownames(marker_file) = marker_file$gene_ID
      
      avg_FC_names = colnames(marker_file)[str_detect(colnames(marker_file), pattern = "avg_log2FC")]
      
      if (nrow(marker_file) == 0) {
        print("No TFs")} 
      
      else if (nrow(marker_file) <= find_candidates & nrow(marker_file) != 0){
        marker_file = marker_file
      } 
      
      else {
        if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) < find_candidates & 
            nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) > 0){
          subset1 = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
          subset2 = marker_file[marker_file[[avg_FC_names[1]]] < 1 & marker_file[[avg_FC_names[2]]] < 1, ]
          subset2 = subset2[order(subset2[, avg_FC_names[1]], subset2[, avg_FC_names[2]], decreasing = TRUE), ]
          
          if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
            marker_file = rbind.data.frame(subset1, subset2)
          } 
          
          else {
            subset2 = head(subset2, (find_candidates - nrow(subset1)))
            marker_file = rbind.data.frame(subset1, subset2)
          }
        }
        else if (nrow(marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]) == 0){
          marker_file = marker_file[order(marker_file[, avg_FC_names[1]], marker_file[, avg_FC_names[2]], decreasing = TRUE), ]
          marker_file = marker_file[c(1:find_candidates), ]
        }
        
        else {
          marker_file = marker_file[marker_file[[avg_FC_names[1]]] >= 1 & marker_file[[avg_FC_names[2]]] >= 1, ]
          marker_file = marker_file[c(1:find_candidates), ]
        }
      }
      
      marker_file$source = "DEGs"
    }
    
    return(marker_file)
  }
  
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
  
  
  # create an empty list to store the candidate markers files
  temp_list = list()
  
  
  for (i in c(1:length(temp_deg_list))) {
    
    # Get cluster specific markers for the cell clusters
    cluster_marker = specific_marker_finder(DEG_file_dir = temp_deg_list[[i]], pct_detection = pct_detection, pct_diff = pct_diff, return_markers_list = TRUE)[[1]]
    
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

}

if (store_outputs == TRUE) {
  
}
