# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("~/Documents/Yasirs_projects/scExplorer/Functions", pattern = "*.R$", full.names = TRUE)
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("~/Documents/Beginning_of_a_compendium/Chapter_2/comparative_study_of_two_species/Liger_Analysis_strategy_2/Seurat_object/seurat_object_of_K_50.RData")

# Cluster ID for which to look for candidate markers
cluster_ID = "15"

# Reduction name storing the basis matrix which will be used to rank the genes of the GEP based on the coefficient
reduction_name = "inmf"

# If GEP ID or IDs are not provided, use the characteristic GEP of the cell cluster to identify candidate markers
defining_GEP = characteristic_GEP_of_cells(seurat_object = integrated.data, target_ID = cluster_ID, store_output = FALSE)

as.character(defining_GEP[which.max(defining_GEP$cell_count), "GEP"])

# The variable can be either a interger or character of GEP ID or a vector of genes
# File path / GEP_ID / vector of genes

# GEP ID from which to extract the candidate markers
GEP_ID = 28

# GEP genes
GEPgenes = read.delim("~/Documents/Yasirs_projects/Analysis_of_single_species_Cardamine/Analysis_with_Liger/Analysis_output/Factor_Ks_GEPs/Factor_K_44_GEPs_Ortho/GEP_GEP_7_ortho.txt", header = FALSE, col.names = "gene_ID")

# Differential expression test to perform. By default, uses the "wilcox" test provided by the Seurat tool

# Run differential expression test for the genes in the desired GEP.
cluster_GEP_markers = FindMarkers(integrated.data, ident.1 = target_ID, features = GEPgenes$gene_ID, test.use = "wilcox", only.pos = TRUE)

# Add the gene names to the DE file
cluster_GEP_markers$gene_ID = rownames(cluster_GEP_markers)

# extract the basis matrix from the Seurat object
basis = as.data.frame(integrated.data@reductions[[reduction_name]]@feature.loadings)

# Add column names to the basis matrix
colnames(basis) = str_c("GEP_", parse_number(colnames(basis)))

basis = basis[GEPgenes, GEP_ID]


# Create a list to store the GEP genes, as this will facilitate the iteration if there is more than one GEPs


# Order the matrix with all geps
basis = order_factorized_matrix(basis)

Basis_subset = basis[GEPgenes$gene_ID, GEP_ID, drop = FALSE]
colnames(Basis_subset)[1] <- "coefficient"
Basis_subset = Basis_subset[order(Basis_subset$coefficient, decreasing = TRUE), , drop = FALSE]
Basis_subset$rank = c(1:nrow(Basis_subset))
Basis_subset$gene_ID= rownames(Basis_subset)

# lets create a description file for the basis genes
description_basis_genes = ATH[rownames(ATH) %in% rownames(Basis_subset), ]

# Description - All DEGs identified for the cluster
geps_description = merge(Basis_subset, description_basis_genes, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(5, 6)]

geps_description = geps_description[order(geps_description$rank, decreasing = FALSE), ]

geps_description = merge(geps_description, cluster_GEP_markers, by = "gene_ID")

geps_description = geps_description[, c(1, 4, 7, 10, 11, 12, 2, 3, 9, 13, 5, 6, 8)]

geps_description = geps_description[order(geps_description$pct.2, decreasing = FALSE), ]

# if the selection criteria does not work, we picked the top 10 genes
geps_description = geps_description[geps_description$pct.2 < 0.1, ]

if (nrow(geps_description) == 0) {
  print("No genes")} else if (nrow(geps_description) <= 10 & nrow(geps_description) != 0){
    geps_description = geps_description
  } else {
    if (nrow(geps_description[geps_description$avg_log2FC >= 1,]) < 10 & nrow(geps_description[geps_description$avg_log2FC >= 1,]) > 0){
      subset1 = geps_description[geps_description$avg_log2FC >= 1,]
      subset2 = geps_description[geps_description$avg_log2FC < 1,]
      subset2 = subset2[order(subset2$avg_log2FC, decreasing = TRUE), ]
      
      if (nrow(subset2) <= (10 - nrow(subset1))) {
        geps_description = rbind.data.frame(subset1, subset2)
      } else {
        subset2 = head(subset2, (10 - nrow(subset1)))
        geps_description = rbind.data.frame(subset1, subset2)
      }
    }
    else if (nrow(geps_description[geps_description$avg_log2FC >= 1,]) == 0){
      geps_description = geps_description[order(geps_description$avg_log2FC, decreasing = TRUE), ]
      geps_description = geps_description[c(1:10), ]
    }
    else {
      geps_description = geps_description[geps_description$avg_log2FC >= 1, ]
      geps_description = geps_description[c(1:10), ]
    }
  }

rownames(geps_description) <- geps_description$gene_ID

geps_description$source = str_c("GEP_", GEP_ID)

write.csv(geps_description, file = str_c("Cluster_", target_ID, "_GEP_", GEP_ID,"_markers.csv"), FALSE)

if (nrow(geps_description) < 10){
  
  subset1 = geps_description[ , !colnames(geps_description) %in% "source"]
  
  geps_description = merge(Basis_subset, description_basis_genes, by.x = "gene_ID", by.y = "CH_ID", all.x = TRUE)[-c(5, 6)]
  
  geps_description = geps_description[order(geps_description$rank, decreasing = FALSE), ]
  
  geps_description = merge(geps_description, cluster_GEP_markers, by = "gene_ID")
  
  geps_description = geps_description[, c(1, 4, 7, 10, 11, 12, 2, 3, 9, 13, 5, 6, 8)]
  
  geps_description = geps_description[order(geps_description$pct.2, decreasing = FALSE), ]
  
  subset2 = geps_description[!geps_description$gene_ID %in% rownames(subset1), ]
  
  if (nrow(subset2) <= (10 - nrow(subset1))) {
    geps_description = rbind.data.frame(subset1, subset2)
  } else {
    subset2 = head(subset2, (10 - nrow(subset1)))
    geps_description = rbind.data.frame(subset1, subset2)
  }
  
}


rownames(geps_description) <- geps_description$gene_ID

geps_description$source = str_c("GEP_", GEP_ID)

write.csv(geps_description, file = str_c("Cluster_", target_ID, "_GEP_", GEP_ID,"_markers.csv"), FALSE)


