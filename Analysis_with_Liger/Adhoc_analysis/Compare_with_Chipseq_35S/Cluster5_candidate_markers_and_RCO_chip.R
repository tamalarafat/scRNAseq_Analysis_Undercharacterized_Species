# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_3/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

if (!dir.exists("Genes_peak")){
  dir.create("Genes_peak", showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Load the chipseq genes
genes_chipseq = read.delim("idr.optimal_peak.narrowPeak_RES50_35SRCO_filtered_anno_Ch-AtID_sort.txt")
genes_chipseq$CH_ID = str_replace(genes_chipseq$Chi_ID, "_", "-")

geneID_chipseq = str_replace(genes_chipseq$Chi_ID, "_", "-")

# Keeping only the unique ids
geneID_chipseq = unique(geneID_chipseq)

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# candidate markers
Cluster_markers = read.csv("Cluster_5_candidate_biomarkers.csv", row.names = 1)[, 1]

length(intersect(geneID_chipseq, Cluster_markers)) # 26 (79% match)


# lets create a description file for the basis genes
Cluster_markers_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster_markers), ]

# Write a csv file with Chipseq
write.csv(Cluster_markers_chip_subset, file = "Genes_peak/Chip_subset_cluster5_candidate_markers.csv", row.names = FALSE)

Cluster_markers_chip_subset = Cluster_markers_chip_subset[Cluster_markers_chip_subset$Annotation == "promoter-TSS", ]

# Write a csv file with Chipseq
write.csv(Cluster_markers_chip_subset, file = "Genes_peak/Cluster_5_Final_Chip_subset.csv", row.names = FALSE)

length(unique(Cluster_markers_chip_subset$CH_ID)) # 23

# Get the description - candidate markers (SM)
description_SM_and_Chip = ATH[rownames(ATH) %in% intersect(geneID_chipseq, Cluster_markers), ]

temp_mat = matrix("NA", ncol = ncol(description_SM_and_Chip), nrow = length(setdiff(intersect(geneID_chipseq, Cluster_markers), rownames(description_SM_and_Chip))))
temp_mat = data.frame(temp_mat)
colnames(temp_mat) = colnames(description_SM_and_Chip)

for (i in c(1: length(setdiff(intersect(geneID_chipseq, Cluster_markers), rownames(description_SM_and_Chip))))){
  temp_mat$CH_ID[i] = setdiff(intersect(geneID_chipseq, Cluster_markers), rownames(description_SM_and_Chip))[i]
}

temp = rbind.data.frame(description_SM_and_Chip, temp_mat)
rownames(temp) = temp$CH_ID

# Write a csv file with Chipseq
write.csv(temp, file = "Genes_peak/Cluster5_candidate_markers_Chip_genes_description.csv", row.names = FALSE)


