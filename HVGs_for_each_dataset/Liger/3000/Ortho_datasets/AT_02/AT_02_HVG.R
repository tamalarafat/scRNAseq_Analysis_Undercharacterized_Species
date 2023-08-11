# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###

###
# Protoplasting induced genes
###

# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = readLines("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/AT_orthologues.txt")

###
# WT C. hirsuta
###

# Load data - WT COL 1st Experiment (leaf 5 and 6)
COL_data <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col_RNA_2nd_ALL_2/filtered_feature_bc_matrix/")

# remove the missing genes from the data
COL_data <- COL_data[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_data)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)

##### Remove the protoplasting induced genes
COL_data <- COL_data[genes_to_keep, ]

###
# COL - 2 E
###

# First replicate - COL 2E - total cells 10760; filter out genes that are not detected in at least 21 cells
COL_2E <- CreateSeuratObject(counts = COL_data, project = "COL_2E", min.features = 200)

# Add metadata information to the seurat object
COL_2E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Thaliana", "WT-COL-2", "WT", "Leaf")

# Remove cells with a total count more than 110000
COL_2E <- subset(COL_2E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
COL_2E[["percent.mt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATM")

# calculate the percentage of total counts belonging to the chloroplast genes.
COL_2E[["percent.pt"]] <- PercentageFeatureSet(COL_2E, pattern = "^ATC")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
COL_2E <- subset(COL_2E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize data to relative count as in Liger
COL_2E <- NormalizeData(COL_2E, normalization.method = "RC", scale.factor = 1)


# Selection of Highly Variable Genes (Liger)
COL_HVGs = variableGenes_liger(dataObject = COL_2E, var.thresh = 0.1, num.genes = 3000, alpha.thresh = 0.99, tol = 0.0001)
fileGenerator(COL_HVGs, fileName = "AT_2E_3000_HVG.txt")

mean_var_COL_2E = generate_mean_var_liger(COL_2E)

save(mean_var_COL_2E, file = "Mean_Variance_AT_02.RData")

# Define color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

p <- ggplot(data = mean_var_COL_2E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_COL_2E[COL_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  geom_abline(intercept = log10(mean((1 / COL_2E@meta.data[["nCount_RNA"]]))), slope = 1, color = "red") + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "COL_2E_Liger_HVG_ab.png", width = 12, height = 12, dpi = 300)

# Without abline
p <- ggplot(data = mean_var_COL_2E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_COL_2E[COL_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "COL_2E_Liger_HVG.png", width = 12, height = 12, dpi = 300)

