# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

HVGs = integrated.data@assays$RNA@var.features

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

hvg_subset = ATH[rownames(ATH) %in% HVGs, ]

# load the cell cycle genes
ccg = read.csv("List_Cell_cycle_genes_from_various_datasets.csv")
ccg = toupper(ccg$ID)

GEP8 = read.delim("GEP_8_ortho_genes.txt", header = FALSE) 

GEP11 = read.delim("GEP_11_ortho_genes.txt", header = FALSE) 

GEP25 = read.delim("GEP_25_ortho_genes.txt", header = FALSE) 

GEP44 = read.delim("GEP_44_ortho_genes.txt", header = FALSE) 

ccg_in_GEPs_ortho_genotypes = c(GEP8$V1, GEP11$V1, GEP25$V1, GEP44$V1)
writeLines(ccg_in_GEPs_ortho_genotypes, "ccg_in_GEPs_ortho_genotypes.txt")

length(intersect(GEP8$V1, ccg)) # = 0

length(intersect(GEP11$V1, ccg)) # = 36

length(intersect(GEP25$V1, ccg)) # = 7

length(intersect(GEP44$V1, ccg)) # = 29

length(intersect(hvg_subset$AT_ID, ccg)) # = 0

writeLines(intersect(hvg_subset$AT_ID, ccg), "ccg_in_hvgs_ortho_genotypes.txt")

# genes

GEP8 = read.delim("GEP_8_genes.txt", header = FALSE) 

GEP11 = read.delim("GEP_11_genes.txt", header = FALSE) 

GEP25 = read.delim("GEP_25_genes.txt", header = FALSE) 

GEP44 = read.delim("GEP_44_genes.txt", header = FALSE) 

ccg_in_GEPs_genotypes = c(GEP8$V1, GEP11$V1, GEP25$V1, GEP44$V1)
writeLines(ccg_in_GEPs_genotypes, "ccg_in_GEPs_genotypes.txt")
