# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

col.genes = rownames(integrated.data)

####################################################################################################

# Load the saved integrated file
CH_TI <- loadRData("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Trajectory/WT_CH_TI.RData")

# Lets assign the cell type name in the data
levels(integrated.data)

# Cell clusters are color coded based on their tyepes
legendFill = hue_pal()(length(levels(integrated.data)))
mypalCC = legendFill[integrated.data$RNA_snn_res.0.3]

# Store  the tsne and umap emdbeddings to the slingshotObject
counts <- as.matrix(integrated.data@assays$RNA@data[integrated.data@assays$RNA@var.features, ])

sce <- SingleCellExperiment(assays = List(counts = counts))

#UMAP
CH_TI@reducedDim <- Embeddings(integrated.data, "umap")

reducedDims(sce) <- SimpleList(UMAP = reducedDims(CH_TI))

# UMAP
png("Figure_7A.png", width = 5600, height = 5600, res = 300)
plot(reducedDims(sce)$UMAP, col = mypalCC, asp = 1, pch = 16, cex = 1, bty = "n", yaxt = "n", xaxt = "n", ann = FALSE)
lines(SlingshotDataSet(CH_TI), type = "lineages", lwd = 6, col = "black", bty = "n")
dev.off()

