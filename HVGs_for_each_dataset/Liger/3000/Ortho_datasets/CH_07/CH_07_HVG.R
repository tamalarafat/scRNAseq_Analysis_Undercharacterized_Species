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
OX_data <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Ox_RNA_7th_ALL_2_Newest/filtered_feature_bc_matrix/")

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
# OX - 7 E
###

# Fourth replicate - OX 7E - total cells 9090; filter out genes that are not detected in at least 18 cells
OX_7E <- CreateSeuratObject(counts = OX_ortho_data, project = "OX_7E", min.features = 200)

# Add metadata information to the seurat object
OX_7E[[c("Species", "Replicates", "Genotype", "Tissue")]] <- c("Hirsuta", "WT-OX-7", "WT", "Leaf")

# Remove cells with a total count more than 110000
OX_7E <- subset(OX_7E, subset = nCount_RNA <= 110000)

# calculate the percentage of total counts belonging to the mitochondiral genes.
OX_7E[["percent.mt"]] <- PercentageFeatureSet(OX_7E, pattern = "Mt")

# calculate the percentage of total counts belonging to the chloroplast genes.
OX_7E[["percent.pt"]] <- PercentageFeatureSet(OX_7E, pattern = "Pt")

# Remove cells using the mitochondiral percentage and chloroplast percentage threshold
OX_7E <- subset(OX_7E, subset = percent.mt < 5 & percent.pt < 10)

# Normalize data to relative count as in Liger
OX_7E <- NormalizeData(OX_7E, normalization.method = "RC", scale.factor = 1)


# Selection of Highly Variable Genes (Liger)
OX_HVGs = variableGenes_liger(dataObject = OX_7E, var.thresh = 0.1, num.genes = 3000, alpha.thresh = 0.99, tol = 0.0001)
fileGenerator(OX_HVGs, fileName = "CH_7E_3000_HVG.txt")

mean_var_OX_7E = generate_mean_var_liger(OX_7E)

save(mean_var_OX_7E, file = "Mean_Variance_CH_07.RData")

# Define color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

p <- ggplot(data = mean_var_OX_7E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_7E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  geom_abline(intercept = log10(mean((1 / OX_7E@meta.data[["nCount_RNA"]]))), slope = 1, color = "red") + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "OX_7E_Liger_HVG_ab.png", width = 12, height = 12, dpi = 300)

# Without abline
p <- ggplot(data = mean_var_OX_7E, aes(x = avgExpression, y = geneVariance)) + 
  geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = mean_var_OX_7E[OX_HVGs, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggsave(plot = p, filename = "OX_7E_Liger_HVG.png", width = 12, height = 12, dpi = 300)

