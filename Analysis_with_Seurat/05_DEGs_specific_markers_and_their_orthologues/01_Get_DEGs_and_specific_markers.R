# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Projects_Yasir/scExplorer/Functions/", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

# Load the gene description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Dir - containing the DEG files
DEG_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/Differentially_expressed_genes/DEGs_of_each_cluster/DEG_DEtest_wilcox"

markers_list = get_DEGs_details(DEG_file_dir = DEG_dir, gene_description_file = ATH, return_markers_list = TRUE)

save(markers_list, file = "DEGs_and_Markers.RData")

