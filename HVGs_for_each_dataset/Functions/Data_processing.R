rowMeansFast <- function(x) {
  .Call('_rliger_rowMeansFast', PACKAGE = 'rliger', x)
}

rowVarsFast <- function(x, means) {
  .Call('_rliger_rowVarsFast', PACKAGE = 'rliger', x, means)
}

SparseRowVar2 <- function(mat, mu, display_progress) {
  .Call('_Seurat_SparseRowVar2', PACKAGE = 'Seurat', mat, mu, display_progress)
}

SparseRowVarStd <- function(mat, mu, sd, vmax, display_progress) {
  .Call('_Seurat_SparseRowVarStd', PACKAGE = 'Seurat', mat, mu, sd, vmax, display_progress)
}

FastExpMean <- function(mat, display_progress) {
  .Call('_Seurat_FastExpMean', PACKAGE = 'Seurat', mat, display_progress)
}

FastLogVMR <- function(mat, display_progress) {
  .Call('_Seurat_FastLogVMR', PACKAGE = 'Seurat', mat, display_progress)
}



# Load the table containing the list of protoplasting-induced genes.
PP_genes_table = read.csv("/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/Protoplasting_genes/Ox_Co0_leaf_protoplast_v12_final_August_2022_ortho.csv")

# Gene IDs - protoplasting-induced genes
PP_genes = PP_genes_table$GeneID

###
# Orthologues table
###

RCO = "AT5G67651"

ortho_table = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Input_files/Additional_inputs/Orthologues_n_correspondence/Orthos_table.csv")

###
# WT C. hirsuta
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
OX_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_OX_RNA_1ST_2_Newest/filtered_feature_bc_matrix/")

# Convert the gene ids in the data table to ortho gene ids
OX_DF_1 = prepare_ortho_data(input_data = OX_data_1E, ortho_data = ortho_table, ortho_column_name_of_gene_ids = "C.hirsutaOX", ortho_column_name_to_assign = "A.thaliana.TAIR10")

###
# WT A. thaliana
###

# Load data - WT OX 1st Experiment (leaf 5 and 6)
COL_data_1E <- Read10X(data.dir = "/netscratch/dep_tsiantis/common/scRNAseq/FINAL_datasets_200822/outs_Col0_RNA_1ST_2/filtered_feature_bc_matrix/")

# All gene IDs - Arabidopsis Thaliana
thaliana_genes = rownames(COL_data_1E)

RCO %in% thaliana_genes

# extracting the Cardamine IDs that are present in orthologues table 
thaliana_ortho_genes = as.character(ortho_table$A.thaliana.TAIR10)

# not all the thaliana ids are present in the ortho data - 51 thaliana genes are missing in the thaliana data
thaliana_ortho_genes = intersect(thaliana_genes, thaliana_ortho_genes)

# Let's subset the data with the ortho genes
COL_data_1E <- COL_data_1E[thaliana_ortho_genes, ]

# remove the missing genes from the data
OX_DF_1 <- OX_DF_1[thaliana_ortho_genes, ]

# All gene IDs - all datasets
ortho_genes = rownames(COL_data_1E)

# Remove protoplasting-induced genes from the total set of hirsuta genes
genes_to_keep = setdiff(ortho_genes, PP_genes)

##### Remove the protoplasting induced genes
OX_DF_1 <- OX_DF_1[genes_to_keep, ]

###
# OX - 1 E
###

# First replicate - OX 1E - total cells 6640; filter out genes that are not detected in at least 13 cells
OX_1E <- CreateSeuratObject(counts = OX_DF_1, project = "OX_1E", min.features = 200)

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

# Normalize data to relative count as in Liger
OX_1E <- NormalizeData(OX_1E, normalization.method = "RC", scale.factor = 1)

# OX_1E <- FindVariableFeatures(OX_1E, selection.method = "vst", nfeatures = 2000)

# Extract the count matrix
OX_W1 <- GetAssayData(OX_1E, assay = "RNA", slot = "counts")

# Add a string (identifier) with the barcodes
colnames(OX_W1) <- paste("O1", colnames(OX_W1), sep = "_")

