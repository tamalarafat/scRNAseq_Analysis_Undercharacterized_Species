# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Seurat_Analyses/Seurat_analysis_strategy_1/1_Integration_script/integrated_ox_wt_seurat.RData")

Idents(integrated.data) <- "integrated_snn_res.0.4"

integrated.data$reps <- integrated.data$orig.ident
integrated.data$reps <- factor(integrated.data$reps, levels = c("OX_1E", "OX_2E", "OX_3E", "OX_7E"), labels = c("O1", "O2", "O3", "O7"))

# Get the metadata file
md = integrated.data@meta.data
md$cell_ID = rownames(md)

md = md[, c("cell_ID", "integrated_snn_res.0.4", "reps")]

md$cell_ID = str_c(md$reps, "_", substr(md$cell_ID, start = 1, stop = 18))
md = md[, c(1, 2)]
colnames(md)[2] <- "seurat_withCC"

save(md, file = "seurat_withCC.RData")
