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
OX_data <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_3rd_ALL_3000_Newest/filtered_feature_bc_matrix/")

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
# OX - 3 E
###

# Third replicate - OX 3E - total cells 4100; filter out genes that are not detected in at least 8 cells
OX_3E <- CreateSeuratObject(counts = OX_ortho_data, project = "OX_3E", min.features = 200)

# Add metadata information to the seurat object
OX_3E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-3", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_3E <- subset(OX_3E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_3E[["percent.mt"]] <- PercentageFeatureSet(OX_3E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_3E[["percent.pt"]] <- PercentageFeatureSet(OX_3E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_3E <- subset(OX_3E, subset = percent.mt < 5 & percent.pt < 10)
# Normalize data to relative count as in Liger
OX_3E <- NormalizeData(OX_3E, normalization.method = "RC", scale.factor = 1)


# Selection of Highly Variable Genes (Liger)
OX_HVGs = variableGenes_liger(dataObject = OX_3E, var.thresh = 0.1, num.genes = 3000, alpha.thresh = 0.99, tol = 0.0001)
fileGenerator(OX_HVGs, fileName = "CH_3E_3000_HVG.txt")

mean_var_OX_3E = generate_mean_var_liger(OX_3E)

save(mean_var_OX_3E, file = "Mean_Variance_CH_03.RData")

# Define color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

p <- ggplot(data = mean_var_OX_3E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_3E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  geom_abline(intercept = log10(mean((1 / OX_3E@meta.data[["nCount_RNA"]]))), slope = 1, color = "red") + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "OX_3E_Liger_HVG_ab.png", width = 12, height = 12, dpi = 300)

# Without abline
p <- ggplot(data = mean_var_OX_3E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_3E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "OX_3E_Liger_HVG.png", width = 12, height = 12, dpi = 300)

