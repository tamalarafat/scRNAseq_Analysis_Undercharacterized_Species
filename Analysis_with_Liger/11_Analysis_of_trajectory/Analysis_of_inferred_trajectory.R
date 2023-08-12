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

## UMAP plot smooth - the smooth curvve does not work
## UMAP plot 1200 * 1200 -  UMAP_Straight_CC_TSNE

#TSNE
CH_TI@reducedDim <- Embeddings(integrated.data, "tsne")

reducedDims(sce) <- SimpleList(TSNE = reducedDims(CH_TI))

# Lets save the plots in folder system
dir.create("TI_PCA", showWarnings = TRUE, recursive = FALSE, mode = "0777")

# Global Lineage 

dir.create("TI_PCA/Reduced_Dimensions", showWarnings = T, recursive = F, mode = "0777")

#UMAP
CH_TI@reducedDim <- Embeddings(integrated.data, "umap")

reducedDims(sce) <- SimpleList(UMAP = reducedDims(CH_TI))

# UMAP
png("TI_PCA/Reduced_Dimensions/UMAP_TI_onPCA_CCL.png", width = 4800, height = 4800, res = 300)
plot(reducedDims(sce)$UMAP, col = mypalCC, asp = 1, pch = 16, cex = 1, bty = "n", yaxt = "n", xaxt = "n", ann = FALSE)
lines(SlingshotDataSet(CH_TI), type = "lineages", lwd = 3, col = "black", bty = "n")
dev.off()

# Lets get the cells in each lineage by plotting the corresponding curve
NC <- 3 # Number of columns
PT <- slingPseudotime(CH_TI) # Pseudotime ordering values of the cells
PT <- cbind(PT, Seurat = as.numeric(integrated.data$RNA_snn_res.0.3))
nms <- colnames(PT) # column names of the pseudotime ordering table
NR <- ceiling(length(nms)/NC) # Number of rows
pal <- viridis(100, end = 0.95) # Palette to use colors for each cell

# Curve 1
summary(CH_TI@curves)

# Let's assign UMAP embeddings to the Slingshot object
#UMAP
CH_TI@reducedDim <- Embeddings(integrated.data, "umap")

reducedDims(sce) <- SimpleList(UMAP = reducedDims(CH_TI))

dir.create("TI_PCA/Lineage_Structures", showWarnings = T, recursive = F, mode = "0777")

# Lets use the TSNE plot for visualization
# Lineage 1
paste(c("Lineage 1 :", CH_TI@lineages$Lineage1), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage1.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 1], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 1 : cell clusters"), paste(CH_TI@lineages$Lineage1, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()


# Lineage 2
paste(c("Lineage 2 :", CH_TI@lineages$Lineage2), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage2.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 2], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 2 : cell clusters"), paste(CH_TI@lineages$Lineage2, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()


# Lineage 3
paste(c("Lineage 3 :", CH_TI@lineages$Lineage3), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage3.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 3], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 3 : cell clusters"), paste(CH_TI@lineages$Lineage3, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()


# Lineage 4
paste(c("Lineage 4 :", CH_TI@lineages$Lineage4), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage4.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 4], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 4 : cell clusters"), paste(CH_TI@lineages$Lineage4, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()

# Lineage 5
paste(c("Lineage 5 :", CH_TI@lineages$Lineage5), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage5.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 5], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 5 : cell clusters"), paste(CH_TI@lineages$Lineage5, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()

# Lineage 6
paste(c("Lineage 6 :", CH_TI@lineages$Lineage6), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage6.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 6], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 6 : cell clusters"), paste(CH_TI@lineages$Lineage6, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()

# Lineage 7
paste(c("Lineage 7 :", CH_TI@lineages$Lineage7), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage7.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 7], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 7 : cell clusters"), paste(CH_TI@lineages$Lineage7, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()

# Lineage 8
paste(c("Lineage 8 :", CH_TI@lineages$Lineage8), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage8.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 8], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 8 : cell clusters"), paste(CH_TI@lineages$Lineage8, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()


# Lineage 9
paste(c("Lineage 9 :", CH_TI@lineages$Lineage9), collapse = "  ")

png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_Lineage9.png", width = 2800, height = 2800, res = 300)
colors <- pal[cut(PT[, 9], breaks = 100)]
plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n",
     main = paste(c("Lineage 9 : cell clusters"), paste(CH_TI@lineages$Lineage9, collapse = " > ")))
lines(SlingshotDataSet(CH_TI), type = "lineages", col = "black", lwd = 2)
dev.off()


# Lets put all the curves in one figure
png("TI_PCA/Lineage_Structures/UMAP_TI_onPCA_AllLineages.png", width = 2800, height = 2800, res = 300)
par(mfrow = c(NR, NC))

for (i in nms) {
  colors <- pal[cut(PT[,i], breaks = 100)]
  plot(reducedDim(CH_TI), col = colors, pch = 16, cex = 0.5, asp = 1, bty = "n", main = i)
  lines(CH_TI, lwd = 2, col = 'black', type = 'lineages')
}

dev.off()

############################################################################################

dir.create("TI_PCA/Pseudotime_Values_ofCells", showWarnings = T, recursive = F, mode = "0777")

### Lin 1
paste(c("Lineage 1 :", CH_TI@lineages$Lineage1), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,1])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "3", "17", "0", "10")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "3", "17", "0", "10"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 1 ordered by pseudotime") + 
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))

ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage1.png", 12, 12)                    


### Lin 2
paste(c("Lineage 2 :", CH_TI@lineages$Lineage2), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,2])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "3", "17", "0", "14")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "3", "17", "0", "14"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 2 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))

ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage2.png", 12, 12)   


### Lin 3
paste(c("Lineage 3 :", CH_TI@lineages$Lineage3), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,3])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "3", "15", "8")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "3", "15", "8"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 3 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage3.png", 12, 12)   


### Lin 4
paste(c("Lineage 4 :", CH_TI@lineages$Lineage4), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,4])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "2", "12")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "2", "12"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 4 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage4.png", 12, 12)   


### Lin 5
paste(c("Lineage 5 :", CH_TI@lineages$Lineage5), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,5])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "4", "1")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "4", "1"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 5 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage5.png", 12, 12)   


### Lin 6
paste(c("Lineage 6 :", CH_TI@lineages$Lineage6), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,6])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "2", "7")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "2", "7"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 6 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage6.png", 12, 12)  


### Lin 7
paste(c("Lineage 7 :", CH_TI@lineages$Lineage7), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,7])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "4", "11")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "4", "11"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 7 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage7.png", 12, 12)  


### Lin 8
paste(c("Lineage 8 :", CH_TI@lineages$Lineage8), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,8])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "5")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "5"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 8 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage8.png", 12, 12)  


### Lin 9
paste(c("Lineage 9 :", CH_TI@lineages$Lineage9), collapse = "  ")
slingshot_df = data.frame(Seurat = integrated.data$RNA_snn_res.0.3, slingPseudotime_1 = PT[,9])
slingshot_df = slingshot_df[!is.na(slingshot_df$slingPseudotime_1), ]
cells_to_Keep = c("13", "6", "9", "16")
slingshot_df = slingshot_df[slingshot_df$Seurat %in% cells_to_Keep, ]
slingshot_df$Seurat = droplevels(slingshot_df$Seurat)

slingshot_df$Seurat = factor(slingshot_df$Seurat, levels = c("13", "6", "9", "16"))
slingshot_df = slingshot_df[!is.na(slingshot_df$Seurat), ]

p <- ggplot(slingshot_df, aes(x = slingPseudotime_1, y = Seurat, colour = Seurat)) + geom_quasirandom(groupOnX = FALSE)
p <- p +  xlab("Slingshot pseudotime") + ylab("Cell clusters") +
  ggtitle("Cells in lineage 9 ordered by pseudotime") + 
  
  theme(axis.text = element_text(face = "bold", colour = "black", size = 24),
        axis.title = element_text(face = "bold", size = 24),
        legend.key.size = unit(2, "cm"),
        legend.key = element_rect(size = 14),
        legend.text = element_text(size = 24, face = "bold"),
        legend.title = element_text(size = 24, face = "bold"),
        title = element_text(size = 14, face = "bold")) + 
  guides(colour = guide_legend(title = "Clusters", override.aes = list(size = 6)))
ggbasicSaver(p, "TI_PCA/Pseudotime_Values_ofCells/PseudoOrder_Lineage9.png", 12, 12)  

