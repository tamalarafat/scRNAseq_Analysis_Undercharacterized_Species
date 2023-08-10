# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

if (!dir.exists("Genes_peak")){
  dir.create("Genes_peak", showWarnings = TRUE, recursive = FALSE, mode = "0777")
}

# Load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Load the chipseq genes
genes_chipseq = read.delim("idr.optimal_peak.narrowPeak_RES50_35SRCO_filtered_anno_Ch-AtID_sort.txt")
genes_chipseq$CH_ID = str_replace(genes_chipseq$Chi_ID, "_", "-")

geneID_chipseq = str_replace(genes_chipseq$Chi_ID, "_", "-")

# Keeping only the unique ids
geneID_chipseq = unique(geneID_chipseq)

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Specific markers
Cluster5_SM = read.delim("Cluster_5_specific_markers_CH_ID.txt", header = FALSE)[, 1]

length(intersect(geneID_chipseq, Cluster5_SM)) # 87 (79% match)


# lets create a description file for the basis genes
cluster5_SM_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster5_SM), ]

# Write a csv file with Chipseq
write.csv(cluster5_SM_chip_subset, file = "Genes_peak/Chip_subset_cluster5_specific_markers.csv", row.names = FALSE)

# Get the description - specific markers (SM)
description_SM_and_Chip = ATH[rownames(ATH) %in% intersect(geneID_chipseq, Cluster5_SM), ]

temp_mat = matrix("NA", ncol = ncol(description_SM_and_Chip), nrow = length(setdiff(intersect(geneID_chipseq, Cluster5_SM), rownames(description_SM_and_Chip))))
temp_mat = data.frame(temp_mat)
colnames(temp_mat) = colnames(description_SM_and_Chip)

for (i in c(1: length(setdiff(intersect(geneID_chipseq, Cluster5_SM), rownames(description_SM_and_Chip))))){
  temp_mat$CH_ID[i] = setdiff(intersect(geneID_chipseq, Cluster5_SM), rownames(description_SM_and_Chip))[i]
}

temp = rbind.data.frame(description_SM_and_Chip, temp_mat)
rownames(temp) = temp$CH_ID

# Write a csv file with Chipseq
write.csv(temp, file = "Genes_peak/Cluster5_specific_markers_Chip_genes_description.csv", row.names = FALSE)


