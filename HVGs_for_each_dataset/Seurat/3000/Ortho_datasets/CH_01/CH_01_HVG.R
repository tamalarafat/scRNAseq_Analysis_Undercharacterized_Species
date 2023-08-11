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

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the data table to ortho gene ids
OX_ortho_data = prepare_ortho_data(input_data = OX_data, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

# remove the missing genes from the data
OX_ortho_data <- OX_ortho_data[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(OX_ortho_data)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)

##### Remove the protoplasting induced genes
OX_ortho_data <- OX_ortho_data[genes_to_keep, ]

###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_ortho_data, project = "OX_1E", min.features = 200)

# Add metadata information to the seurat object
OX_1E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-1", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_1E <- subset(OX_1E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_1E[["percent.mt"]] <- PercentageFeatureSet(OX_1E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_1E[["percent.pt"]] <- PercentageFeatureSet(OX_1E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_1E <- subset(OX_1E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize the data - log-normalization
OX_1E <- NormalizeData(OX_1E, verbose = FALSE)

# Find a set of highly avariable genes - 3000 HVGs
OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 3000)


# Selection of Highly Variable Genes (Seurat)
OX_HVGs = variableGenes_seurat(dataObject = OX_1E, nfeatures = 3000)
fileGenerator(OX_HVGs, fileName = "CH_1E_3000_HVG_seurat.txt")

mean_var_OX_1E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/CH_01/Mean_Variance_CH_01.RData")

seurat_mean_var_OX_1E = get_variableGenes_info(OX_1E)
save(seurat_mean_var_OX_1E, file = "Seurat_Mean_Variance_CH_01.RData")

VariableFeaturePlot(OX_1E)

# Define color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

p <- ggplot(data = mean_var_OX_1E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_1E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  geom_smooth(method = "loess", se = F, aes(group=NULL), color = "green") +
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))


ggsave(plot = p, filename = "OX_1E_Seurat_HVG_sm.png", width = 12, height = 12, dpi = 300)

# Without abline
p <- ggplot(data = mean_var_OX_1E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_1E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "OX_1E_Seurat_HVG.png", width = 12, height = 12, dpi = 300)


