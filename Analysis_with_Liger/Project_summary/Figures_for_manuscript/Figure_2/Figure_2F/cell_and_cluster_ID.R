# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

md = integrated.data@meta.data
md$cell_ID = rownames(md)

md = md[, c("cell_ID", "RNA_snn_res.0.3")]
colnames(md)[2] <- "withAll"

save(md, file = "withAll.RData")

cl14 = md[md$withAll == 14, ]

without = loadRData("/home/yasir/Documents/Thesis_PhD/Chapter_1/Remove_rep_specific_GEP20_and_22/GEPs_distribution/withCC.RData")

cl14cells_in_without = without[cl14$cell_ID, ]

cl19 = md[md$withAll == 19, ]

cl19cells_in_without = without[cl19$cell_ID, ]
