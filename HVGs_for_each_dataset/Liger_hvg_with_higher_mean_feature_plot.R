# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

load("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Analyses/Species_WT/Liger_06_03_2023/Seurat_object_K_50/On_coefficient/seurat_object_of_K_50.RData")

gene_list = read_lines("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Compare_liger_and_seurat/WT_Species_comparison/Liger_hvg_with_higher_mean.txt")

# Storing directory
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Analyses/Species_WT/Liger_06_03_2023/Seurat_object_K_50/On_coefficient"

for (i in c(1:length(gene_list))){
  feature_expression_split_visualization(seuratObject = integrated.data, store_dir = storing_dir, split_variable = "Species", gene_ID = gene_list[i], gene_name = gene_list[i])
}
