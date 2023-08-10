# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# List all the files
marker_files = str_sort(list.files(pattern = "Cluster"), numeric = TRUE)

# Store them in a different folder
if(!dir.exists("Cluster_candidate_markers_final_list")){
  dir.create("Cluster_candidate_markers_final_list", showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

temp_dir = str_c("Cluster_candidate_markers_final_list", "/")

for (i in c(1:length(marker_files))){
  marker_file = read.csv(marker_files[i], row.names = 1)
  
  criteria_1 = marker_file$pct.2 <= 0.1 # pct.2 less than 10%
  
  criteria_2 = (marker_file$pct.1 - marker_file$pct.2) >= 0.3 # Difference between pct1 and pct2 ust be 30% or more cells expressing the gene

  marker_file = marker_file[criteria_1 | criteria_2, ] # if a gene meet one of these two criteria, keep it in the candidate list
  
  marker_file = marker_file[order(marker_file$pct.2), ]
  
  marker_file$rank = c(1:nrow(marker_file)) # ranking based on pct
  
  write.csv(marker_file, file = str_c(temp_dir, "Cluster_", parse_number(marker_files[i]), "_candidate_biomarkers.csv"))
}

