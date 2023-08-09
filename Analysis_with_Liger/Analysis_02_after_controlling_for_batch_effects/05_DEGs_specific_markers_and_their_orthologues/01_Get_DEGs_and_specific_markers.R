# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Chapter1/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Load the seurat object
load("/home/yasir/Documents/Chapter1/Remove_rep_specific_GEP20_and_22/seurat_object_of_K_44_Without_GEP20.RData")

# Dir - containing the DEG files
DEG_dir = "/home/yasir/Documents/Chapter1/Remove_rep_specific_GEP20_and_22/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

markers_list = get_DEGs_details(DEG_file_dir = DEG_dir, gene_description_file = ATH, return_markers_list = TRUE)

save(markers_list, file = "DEGs_and_Markers.RData")
