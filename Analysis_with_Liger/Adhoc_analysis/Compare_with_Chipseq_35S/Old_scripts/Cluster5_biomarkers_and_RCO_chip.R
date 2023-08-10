# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
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

# Biomarkers
Cluster5_biomarkers = read.csv("Cluster_5_candidate_biomarkers.csv", row.names = 1)

length(intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID)) # 26 genes

# lets create a description file for the basis genes
cluster5_biomarkers_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), ]

# Write a csv file with Chipseq
write.csv(cluster5_biomarkers_chip_subset, file = "Genes_peak/Chip_subset_cluster5_candidate_biomarkers.csv", row.names = FALSE)

# Get the description - specific markers (SM)
description_biomarkers_and_Chip = ATH[rownames(ATH) %in% intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), ]

temp_mat = matrix("NA", ncol = ncol(description_biomarkers_and_Chip), nrow = length(setdiff(intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), rownames(description_biomarkers_and_Chip))))
temp_mat = data.frame(temp_mat)
colnames(temp_mat) = colnames(description_biomarkers_and_Chip)

for (i in c(1: length(setdiff(intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), rownames(description_biomarkers_and_Chip))))){
  temp_mat$CH_ID[i] = setdiff(intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), rownames(description_biomarkers_and_Chip))[i]
}

temp = rbind.data.frame(description_biomarkers_and_Chip, temp_mat)
rownames(temp) = temp$CH_ID

# Write a csv file with Chipseq
write.csv(temp, file = "Genes_peak/Cluster5_candidate_biomarkers_Chip_genes_description.csv", row.names = FALSE)


