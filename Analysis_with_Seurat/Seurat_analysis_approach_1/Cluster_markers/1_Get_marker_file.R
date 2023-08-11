# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Dir - containing the DEG files
DEG_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

markers_list = get_DEGs_details(DEG_file_dir = DEG_dir, gene_description_file = ATH, return_markers_list = TRUE)

save(markers_list, file = "DEGs_and_Markers.RData")


# get the file names
DEG_files = str_sort(list.files(path = DEG_dir, pattern = "Cluster"), numeric = TRUE)

# load the DEG file
Cluster_DEGs = loadRData(str_c(DEG_dir, "/", DEG_files[8]))

# Prepare the file of cluster differentially expressed gene
Cluster_DEGs = Cluster_DEGs[order(Cluster_DEGs$pct.2, decreasing = FALSE), ]
Cluster_DEGs$gene_ID = rownames(Cluster_DEGs)

# Here create a function for specific marker selection
Cluster_specific_DEGs = Cluster_DEGs[(Cluster_DEGs$pct.2 <= 0.1) | ((Cluster_DEGs$pct.1 - Cluster_DEGs$pct.2) >= 0.3), ]
Cluster_specific_DEGs = Cluster_specific_DEGs[order(Cluster_specific_DEGs$pct.2, decreasing = FALSE), ]
Cluster_specific_DEGs$gene_ID = rownames(Cluster_specific_DEGs)

# Save the file
write.csv(Cluster_specific_DEGs, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers.csv"), row.names = FALSE)

# Write a text file containing Cardamine gene IDs
write_lines(rownames(Cluster_specific_DEGs), file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_CH_ID.txt"))

CS_DEGs_description = merge(Cluster_specific_DEGs, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(8, 9)]
CS_DEGs_description = CS_DEGs_description[order(CS_DEGs_description$pct.2, decreasing = FALSE), ]
rownames(CS_DEGs_description) = CS_DEGs_description$gene_ID

# Save the file
write.csv(CS_DEGs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_description.csv"), row.names = FALSE)

# Write the orthologues ID file
write_lines(CS_DEGs_description$AT_ID[!is.na(CS_DEGs_description$AT_ID)], file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_markers_orthoID.txt"))

# Description - Cluster specific TFs 
CS_TFs_description = CS_DEGs_description[!is.na(CS_DEGs_description$TF), ]

# Save the file
write.csv(CS_TFs_description, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_description.csv"), row.names = FALSE)

# Write the TF hirsuta ID file
write_lines(CS_TFs_description$gene_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_CH_ID.txt"))

# Write the TF ortho ID file
write_lines(CS_TFs_description$AT_ID, file = str_c(temp_store_dir, "Cluster_", parse_number(DEG_files[i]), "_specific_TF_orthoID.txt"))

CH_temp_cluster_DEG_list[[i]] <- rownames(Cluster_DEGs)
names(CH_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))

CH_specific_cluster_marker_list[[i]] <- rownames(Cluster_specific_DEGs)
names(CH_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))

AT_temp_cluster_DEG_list[[i]] <- All_DEGs_description$AT_ID[!is.na(All_DEGs_description$AT_ID)]
names(AT_temp_cluster_DEG_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))

AT_specific_cluster_marker_list[[i]] <- CS_DEGs_description$AT_ID[!is.na(CS_DEGs_description$AT_ID)]
names(AT_specific_cluster_marker_list)[i] <- str_c("Cluster_", parse_number(DEG_files[i]))
