# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Load the seurat object
load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects/seurat_object_of_K_44.RData")

integrated.data$initial.cluster <- integrated.data$RNA_snn_res.0.3

integrated.data@meta.data <- integrated.data@meta.data[, -c(grep(pattern = "RNA_snn_res", colnames(integrated.data@meta.data)))]

integrated.data@reductions[["inmfcc"]] <-integrated.data@reductions[["inmf"]]
integrated.data@reductions[["inmfcc"]]@cell.embeddings <- integrated.data@reductions[["inmfcc"]]@cell.embeddings[, -c(8, 11, 20, 22, 25, 44)]

integrated.data <- RunUMAP(integrated.data, reduction = "inmfcc", dims = 1:dim(integrated.data@reductions[["inmfcc"]])[2], n.components = 2)
integrated.data <- RunTSNE(integrated.data, reduction = "inmfcc", dims = 1:dim(integrated.data@reductions[["inmfcc"]])[2], check_duplicates = FALSE, dim.embed = 2)

integrated.data <- FindNeighbors(integrated.data, reduction = "inmfcc", dims = 1:dim(integrated.data@reductions[["inmfcc"]])[2])

integrated.data <- FindClusters(integrated.data, resolution = seq(0.1, 1.2, 0.1), n.start = 50, n.iter = 50)

Idents(integrated.data) <- "RNA_snn_res.0.3"

save(integrated.data, file = "seurat_object_of_K_44_Without_GEP20_22_and_CC.RData")


