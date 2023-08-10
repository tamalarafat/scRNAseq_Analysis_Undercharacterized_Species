# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Load the GEP 32 genes
GEPgenes = read.delim("GEP_23_genes.txt", header = FALSE, col.names = "gene_ID")

avg_df = AverageExpression(object = integrated.data, assays = "RNA", slot = "counts")
avg_df = avg_df$RNA

gep_df = as.data.frame(avg_df[GEPgenes$gene_ID, ])
gepAvg_df = data.frame(Coefficient = colMeans(gep_df))
gepAvg_df$Clusters = rownames(gepAvg_df)
gepAvg_df$Genes = "avgEx"
gepAvg_df = gepAvg_df[, c("Genes", "Clusters", "Coefficient")]
gep_df$Genes = rownames(gep_df)

long_mat = gep_df %>% melt(
  id.vars = "Genes",
  variable.name = "Clusters",
  value.name = "Coefficient"
)

df = rbind(long_mat, gepAvg_df)

p <- ggplot(data = df[df$Genes != "avgEx", ], aes(x = Clusters, y = Coefficient, group = Genes)) + 
  geom_line(color = "grey", size = 2) + 
  geom_line(data = df[df$Genes == "avgEx", ], aes(x = Clusters, y = Coefficient), color = grp_col[3], size = 4) + 
  xlab("Cell clusters") + 
  ylab("Average expression") + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    axis.title = element_text(size = 42, colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 42, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 42, colour = "black"))

ggsave(filename = "Figure_3H.png", plot = p, width = 12, height = 12, dpi = 300)



# Let's look at the genes that are expressed outside cluster 5

# Load the gene description file

# Load the description file
ATH = read.csv("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# GEP genes description
description_basis_genes = ATH[rownames(ATH) %in% GEPgenes$gene_ID, ]

# Max cluster - 5 genes : Chir05Ox-b20300.2 (PIP2D), Chir05Ox-b12990.2 (ARR5), Chir04Ox-b16010.2 (AtDRM2)
# Max cluster - 10 genes : Chir05Ox-b20300.2 (PIP2D)
# Max cluster - 16 genes : Chir06Ox-b11590.2 (LSH1), Chir04Ox-b13260.2 (FRK2)
# Max cluster - 17 genes : Chir05Ox-b20300.2 (PIP2D)

p <- ggplot(data = df[df$Genes != "avgEx", ], aes(x = Clusters, y = Coefficient, group = Genes)) + 
  geom_line(color = "grey", size = 2) + 
  geom_line(data = df[df$Genes == "Chir05Ox-b20300.2", ], aes(x = Clusters, y = Coefficient), color = grp_col[7], size = 2) + 
  geom_line(data = df[df$Genes == "Chir05Ox-b12990.2", ], aes(x = Clusters, y = Coefficient), color = grp_col[7], size = 2) + 
  geom_line(data = df[df$Genes == "Chir04Ox-b16010.2", ], aes(x = Clusters, y = Coefficient), color = grp_col[7], size = 2) + 
  geom_line(data = df[df$Genes == "Chir06Ox-b11590.2", ], aes(x = Clusters, y = Coefficient), color = grp_col[4], size = 2) + 
  geom_line(data = df[df$Genes == "Chir04Ox-b13260.2", ], aes(x = Clusters, y = Coefficient), color = grp_col[4], size = 2) + 
  geom_line(data = df[df$Genes == "avgEx", ], aes(x = Clusters, y = Coefficient), color = grp_col[3], size = 4) + 
  xlab("Cell clusters") + 
  ylab("Average expression") + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    axis.title = element_text(size = 42, colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 42, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 42, colour = "black"))

ggsave(filename = "Figure_3H_1.png", plot = p, width = 12, height = 12, dpi = 300)

