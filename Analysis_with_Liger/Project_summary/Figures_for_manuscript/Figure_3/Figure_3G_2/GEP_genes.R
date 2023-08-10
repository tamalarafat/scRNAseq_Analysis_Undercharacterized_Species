# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# The GEP we are interested in
GEP_ID = 23

# Work on the basis
basis = as.data.frame(integrated.data@reductions[["inmfcc"]]@feature.loadings)
colnames(basis) = str_c("GEP_", parse_number(colnames(basis)))

# Order the matrix with all geps
basis = order_factorized_matrix(basis)

Basis_subset = basis[, GEP_ID, drop = FALSE]
Basis_subset = Basis_subset[order(Basis_subset$GEP_23, decreasing = FALSE), , drop = FALSE]
Basis_subset$index = c(1:nrow(Basis_subset))

Basis_label = Basis_subset

# Let's load the GEP 23 Cardamine gene IDs which will be used to subset the Basis label
GEPgenes = read.delim("GEP_23_genes.txt", header = FALSE, col.names = "gene_ID")

Basis_label = Basis_label[GEPgenes$gene_ID, ]
Basis_label = Basis_label[order(Basis_label$GEP_23, decreasing = FALSE), ]
Basis_label$gene_ID = rownames(Basis_label)

# Let's load the description file of the GEP genes
gep_des = read.csv("Cluster_13_GEP_23_markers.csv", row.names = 1)
rownames(gep_des) <- gep_des$CH_ID

# Chir04Ox-b12980.2, Chir05Ox-b12990.2, Chir06Ox-b11590.2, Chir06Ox-b35090.2, Chir06Ox-b37780.2, Chir08Ox-b16470.2, Chir08Ox-b28070.2 - All the LSH gene
Basis_label = merge(Basis_label, gep_des, by.x = "gene_ID", by.y = "gene_ID")
Basis_label = Basis_label[order(Basis_label$coefficient, decreasing = FALSE), ]
rownames(Basis_label) <- Basis_label$gene_ID

p <- ggplot(Basis_subset, aes(x = index, y = GEP_23)) + geom_point(color = grp_col[3], size = 4) +
  ggrepel::geom_text_repel(data = Basis_label, aes(label = Gene_code_name), nudge_x = -800, size = 16, force = 10) +
  xlab("GEP 23 genes") + ylab("Coefficient") +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, colour = "black")
  )
ggsave(filename = "Figure_3G_2.png", plot = p, width = 10, height = 16, dpi = 300)

