# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2023/PhD_projects_Yasir/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# # Load the seurat object
# load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Select_UMAP/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

# Path to the DEG files - differentially genes usign the "Findmarkers" function (Seurat) per cluster
DEG_dir = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

cluster_ID = c(15, 16)

find_candidates = 10

specific_markers_list = specific_marker_finder(DEG_file_dir = DEG_dir, pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE)


# Check if any cluster ID is provided or not. If provided, use the specified ID. If not, get all the cluster IDs and go over each ID.
if (missing(cluster_ID)){
  cluster_ID = levels(Idents(integrated.data))
} else {
  cluster_ID = cluster_ID
}

specific_markers_list = specific_markers_list[sapply(cluster_ID, function(x) {grep(pattern = str_c("^", x, "$", sep = ""), parse_number(names(specific_markers_list)))}, simplify = "vector")]

# If it's a DEG file, get the specific marker for that file

# Inside the loop now

# What ever the DEG file is put them in a list and go over the list item

# Narrow the list of the differentially expressed genes by finding the cluster specific markers

Cluster_marker = specific_marker_finder(DEG_file_dir = DEG_dir, file_name_pattern = "Cluster_", pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE, cluster_ID = 15)

Cluster_marker = specific_marker_finder(DEG_file_dir = DEG_dir, file_name_pattern = "Cluster_", pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE, cluster_ID = cluster_ID)


Cluster_marker = Cluster_marker[[1]]

# Let's define the candidates selection criteria
Cluster_marker = Cluster_marker[order(Cluster_marker$pct.2, decreasing = FALSE), ]

if (nrow(Cluster_marker) == 0) {
  print("No TFs")} else if (nrow(Cluster_marker) <= find_candidates & nrow(Cluster_marker) != 0){
    Cluster_marker = Cluster_marker
  } else {
    if (nrow(Cluster_marker[Cluster_marker$avg_log2FC >= 1,]) < find_candidates & nrow(Cluster_marker[Cluster_marker$avg_log2FC >= 1,]) > 0){
      subset1 = Cluster_marker[Cluster_marker$avg_log2FC >= 1,]
      subset2 = Cluster_marker[Cluster_marker$avg_log2FC < 1,]
      subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
      
      if (nrow(subset2) <= (find_candidates - nrow(subset1))) {
        Cluster_marker = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (find_candidates - nrow(subset1)))
        Cluster_marker = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(Cluster_marker[Cluster_marker$avg_log2FC >= 1,]) == 0){
      Cluster_marker = Cluster_marker[order(Cluster_marker$avg_log2FC, decreasing = TRUE), ]
      Cluster_marker = Cluster_marker[c(1:find_candidates), ]
    }
    else {
      Cluster_marker = Cluster_marker[Cluster_marker$avg_log2FC >= 1, ]
      Cluster_marker = Cluster_marker[c(1:find_candidates), ]
    }
  }


Cluster_marker$source = "DEGs"

write.csv(Cluster_marker, str_c("Cluster_", target_ID, "_top_10_TFs.csv"), row.names = FALSE)




















# This does not work - if user provides the file location.
# Needs addition.
cluster_15 = "/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox/Cluster_15_DEGs.RData"

SM_CL15 = specific_marker_finder(DEG_file_dir = cluster_15, pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE)


# This works - if user provides the loaded DEG file. But the name of the item in the list is defined by the function.
# Needs modification.
cluster_15 = loadRData("/home/ytamal2/Documents/2023/PhD_projects_Yasir/Analysis_of_single_species_Cardamine/Analysis_with_Liger/08_Analysis_02_after_controlling_for_batch_effects/Analysis_outputs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox/Cluster_15_DEGs.RData")

SM_CL15 = specific_marker_finder(DEG_file_dir = cluster_15, pct_detection = 0.1, pct_diff = 0.3, return_markers_list = TRUE)


# If user wants a single element of the list - can get
# Needs addition. Simply add an argument and extract from the list

# For TFs
markers_description_TF_subset = markers_description[!is.na(markers_description$TF), ]

# Particular genes set of interest -  could be TFs or other categories
specify_gene_ID = NULL

# Name of the column to be added into the table. If user provides a list of genes, this column will indicate the category -Yes or No 
incorporate_column_name = NULL

# DEG directory or DEG file

# Cluster IDs to look for candidate markers

# Include GEPs for the clusters. Can take the GEPs file already created or GEP IDs for the cluster. If more than one cluster IDs are provided I can use a list item for the clusters

# Incorporate TFs information. Look for candidate markers for the TFs and for the DEGs

# combine with GEPs if asked

# Incorporate this information
# Find markers only from the given set of genes or from the given full set of differentially expressed genes or for both groups



# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Chapter1/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Load the markers from GEPs
GEP_28_markers = read.csv("/home/yasir/Documents/Chapter1/Liger_analysis_strategy_1/Results/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_15/Top_genes_GEP_identity/GEP_28/Cluster_15_GEP_28_markers.csv", row.names = 1)

GEP_42_markers = read.csv("/home/yasir/Documents/Chapter1/Liger_analysis_strategy_1/Results/Analysis_of_WT_Cardamine/Remove_rep_specific_GEP20_and_22/Cluster_biomarkers/Cluster_15/Top_genes_GEP_identity/GEP_42/Cluster_15_GEP_42_markers.csv", row.names = 1)

GEP_markers = rbind(GEP_28_markers, GEP_42_markers)

GEP_markers = GEP_markers[order(GEP_markers$pct.2, decreasing = FALSE), ]

geps_description = GEP_markers[, -c(7, 8)]

cluster_ids = levels(Idents(integrated.data))

target_ID = "15"



# Create ortho genes list and gene description file
description = ATH[rownames(ATH) %in% rownames(Cluster_marker), ]

# Description - All DEGs identified for the cluster
markers_description = merge(Cluster_marker, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(8, 9)]



if (length(intersect(rownames(geps_description), rownames(markers_description_TF_subset))) != 0){
  
  # For TFs
  markers_description_TF_subset = markers_description[!is.na(markers_description$TF), ]
  
  shared_gene = intersect(rownames(geps_description), markers_description_TF_subset$gene_ID)
  
  geps_description[shared_gene, "source"] = str_c(geps_description[shared_gene, "source"], "and DEGs")
  
  markers_description_TF_subset = markers_description_TF_subset[!markers_description_TF_subset$gene_ID %in% shared_gene, ]
  
  if (nrow(markers_description_TF_subset) == 0) {
    print("No TFs")} else if (nrow(markers_description_TF_subset) <= 10 & nrow(markers_description_TF_subset) != 0){
      markers_description_TF_subset = markers_description_TF_subset
    } else {
      if (nrow(markers_description_TF_subset[markers_description_TF_subset$avg_log2FC >= 1,]) < 10 & nrow(markers_description_TF_subset[markers_description_TF_subset$avg_log2FC >= 1,]) > 0){
        subset1 = markers_description_TF_subset[markers_description_TF_subset$avg_log2FC >= 1,]
        subset2 = markers_description_TF_subset[markers_description_TF_subset$avg_log2FC < 1,]
        subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
        
        if (nrow(subset2) <= (10 - nrow(subset1))) {
          markers_description_TF_subset = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (10 - nrow(subset1)))
          markers_description_TF_subset = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(markers_description_TF_subset[markers_description_TF_subset$avg_log2FC >= 1,]) == 0){
        markers_description_TF_subset = markers_description_TF_subset[order(markers_description_TF_subset$avg_log2FC, decreasing = TRUE), ]
        markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
      }
      else {
        markers_description_TF_subset = markers_description_TF_subset[markers_description_TF_subset$avg_log2FC >= 1, ]
        markers_description_TF_subset = markers_description_TF_subset[c(1:10), ]
      }
    }
  
  rownames(markers_description_TF_subset) <- markers_description_TF_subset$gene_ID
  
  markers_description_TF_subset = markers_description_TF_subset[ , c(1, 7, 10, 3, 4, 5, 2, 6, 8, 9, 11)]
  
  markers_description_TF_subset$source = "DEGs"
  
  write.csv(markers_description_TF_subset, str_c("Cluster_", target_ID, "_top_10_TFs.csv"), row.names = FALSE)
  
}

# For genes 
markers_description_genes_subset = markers_description[is.na(markers_description$TF), ]

if (nrow(markers_description_genes_subset) == 0) {
  print("No genes")} else if (nrow(markers_description_genes_subset) <= 10 & nrow(markers_description_genes_subset) != 0){
    markers_description_genes_subset = markers_description_genes_subset
  } else {
    if (nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) < 10 & nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) > 0){
      subset1 = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]
      subset2 = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC < 1,]
      subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
      
      if (nrow(subset2) <= (10 - nrow(subset1))) {
        markers_description_genes_subset = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (10 - nrow(subset1)))
        markers_description_genes_subset = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) == 0){
      markers_description_genes_subset = markers_description_genes_subset[order(markers_description_genes_subset$avg_log2FC, decreasing = TRUE), ]
      markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
    }
    else {
      markers_description_genes_subset = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1, ]
      markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
    }
  }

rownames(markers_description_genes_subset) <- markers_description_genes_subset$gene_ID

markers_description_genes_subset = markers_description_genes_subset[ , c(1, 7, 10, 3, 4, 5, 2, 6, 8, 9, 11)]

markers_description_genes_subset$source = "DEGs"

write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes.csv"), row.names = FALSE)


if (length(intersect(rownames(geps_description), rownames(markers_description_genes_subset))) != 0){
  
  markers_description_genes_subset = markers_description[is.na(markers_description$TF), ]
  
  shared_gene = intersect(rownames(geps_description), markers_description_genes_subset$gene_ID)
  
  geps_description[shared_gene, "source"] = str_c(geps_description[shared_gene, "source"], "and DEGs")
  
  markers_description_genes_subset = markers_description_genes_subset[!markers_description_genes_subset$gene_ID %in% shared_gene, ]
  
  if (nrow(markers_description_genes_subset) == 0) {
    print("No genes")} else if (nrow(markers_description_genes_subset) <= 10 & nrow(markers_description_genes_subset) != 0){
      markers_description_genes_subset = markers_description_genes_subset
    } else {
      if (nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) < 10 & nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) > 0){
        subset1 = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]
        subset2 = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC < 1,]
        subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
        
        if (nrow(subset2) <= (10 - nrow(subset1))) {
          markers_description_genes_subset = rbind.data.frame(subset1, subset2)
        } else {
          subset2 = head(subset2, (10 - nrow(subset1)))
          markers_description_genes_subset = rbind.data.frame(subset1, subset2)
        }
      }
      else if (nrow(markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1,]) == 0){
        markers_description_genes_subset = markers_description_genes_subset[order(markers_description_genes_subset$avg_log2FC, decreasing = TRUE), ]
        markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
      }
      else {
        markers_description_genes_subset = markers_description_genes_subset[markers_description_genes_subset$avg_log2FC >= 1, ]
        markers_description_genes_subset = markers_description_genes_subset[c(1:10), ]
      }
    }
  
  rownames(markers_description_genes_subset) <- markers_description_genes_subset$gene_ID
  
  markers_description_genes_subset = markers_description_genes_subset[ , c(1, 7, 10, 3, 4, 5, 2, 6, 8, 9, 11)]
  
  markers_description_genes_subset$source = "DEGs"
  
  write.csv(markers_description_genes_subset, str_c("Cluster_", target_ID, "_top_10_genes.csv"), row.names = FALSE)
  
}

cluster_biomarker = rbind(geps_description, markers_description_TF_subset, markers_description_genes_subset)

write.csv(cluster_biomarker, file = str_c("Cluster_", target_ID, "_bio_markers.csv"))

