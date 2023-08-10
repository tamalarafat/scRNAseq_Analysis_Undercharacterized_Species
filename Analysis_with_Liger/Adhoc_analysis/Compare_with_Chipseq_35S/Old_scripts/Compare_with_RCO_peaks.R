# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the chipseq genes
genes_chipseq = read.delim("idr.optimal_peak.narrowPeak_RES50_35SRCO_filtered_anno_Ch-AtID_sort.txt")
genes_chipseq$CH_ID = str_replace(genes_chipseq$Chi_ID, "_", "-")

geneID_chipseq = str_replace(genes_chipseq$Chi_ID, "_", "-")

# Keeping only the unique ids
geneID_chipseq = unique(geneID_chipseq)

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# GEP32 genes
GEP32_genes = read.delim("GEP_32_genes.txt", header = FALSE)[, 1]

length(intersect(geneID_chipseq, GEP32_genes)) # 23 genes

GEP32_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, GEP32_genes), ]

# GEP40 genes
GEP40_genes = read.delim("GEP_40_genes.txt", header = FALSE)[, 1]

length(intersect(geneID_chipseq, GEP40_genes))

GEP40_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, GEP40_genes), ]

# Cluster 5 DEGs
Cluster5_DEGs = read.delim("Cluster_5_all_DEGs_CH_ID.txt", header = FALSE)[, 1]

length(intersect(geneID_chipseq, Cluster5_DEGs))

cluster5_DEGs_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster5_DEGs), ]

# Specific markers
Cluster5_SM = read.delim("Cluster_5_specific_markers_CH_ID.txt", header = FALSE)[, 1]

length(intersect(geneID_chipseq, Cluster5_SM))

cluster5_SM_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster5_SM), ]

# Specific markers
Cluster5_biomarkers = read.csv("Cluster_5_bio_markers.csv", row.names = 1)

length(intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID))

cluster5_biomarkers_chip_subset = genes_chipseq[genes_chipseq$CH_ID %in% intersect(geneID_chipseq, Cluster5_biomarkers$gene_ID), ]

