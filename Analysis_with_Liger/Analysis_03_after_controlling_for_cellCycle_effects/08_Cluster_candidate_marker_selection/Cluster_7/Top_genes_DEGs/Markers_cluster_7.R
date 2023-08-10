# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Load the markers from GEPs
GEP_markers = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Remove_GEP20_22_and_CCGEPs/Cluster_DEGs_Markers_Candidates/Cluster_candidates/Cluster_7/Top_genes_GEP_identity/GEP_13/Cluster_7_GEP_13_markers.csv", row.names = 1)

GEP_markers = GEP_markers[order(GEP_markers$pct.2, decreasing = FALSE), ]

geps_description = GEP_markers[, -c(7, 8)]

cluster_ids = levels(Idents(integrated.data))

target_ID = "7"

# Markers from cluster DEGs
# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Thesis_PhD/Chapter_1/Remove_GEP20_22_and_CCGEPs/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

specific_markers_list = specific_marker_finder(DEG_file_dir = DEG_dir, pct_detection = 0.1, pct_diff = 0.3)

# Cluster 7 marker

names(specific_markers_list)[8]

Cluster_marker = specific_markers_list[[8]]

# Create ortho genes list and gene description file
description = ATH[rownames(ATH) %in% rownames(Cluster_marker), ]

# Description - All DEGs identified for the cluster
markers_description = merge(Cluster_marker, description, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(8, 9)]

markers_description = markers_description[order(markers_description$pct.2, decreasing = FALSE), ]

# For TFs
markers_description_TF_subset = markers_description[!is.na(markers_description$TF), ]

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

