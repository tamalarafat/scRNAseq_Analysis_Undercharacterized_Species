# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF")

load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects/seurat_object_of_K_44.RData")

GEP = read.delim("/home/ytamal2/Documents/2023/Beginning_of_a_compendium/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/Figure_1/Figure_1G/GEP_23_ortho_genes.txt", header = FALSE, col.names = "gene_ID")

ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$AT_ID


# Add gene information
GEP_genes_description = ATH[GEP$gene_ID, -c(3, 4)]

write.csv(GEP_genes_description, file ="GEP23_genes_description.csv")
