# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# information retrieval ::  TAIR / Ensembl / NCBI Entrez gene ids for the marker ids
attributes_to_retrieve = c("tair_symbol", "uniprotswissprot", "entrezgene_id")

Thaliana_genes =  read.delim("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/ATgenes.txt", header = FALSE, col.names = "genes")
species_gene_set <- as.character(Thaliana_genes$genes)

# Factor 10 GEP 10 genes

GEP44 =  read.delim("GEP_44_ortho_genes.txt", header = FALSE, col.names = "genes")
GEP44 <- as.character(GEP44$genes)

# For all the genes in the species of interest, query the Ensembl API, and retrieve ids from biomartr 
species_genes_annotated <- biomartr::biomart(genes = species_gene_set,
                                             mart       = "plants_mart",                 
                                             dataset    = "athaliana_eg_gene",           
                                             attributes = attributes_to_retrieve,        
                                             filters    =  "ensembl_gene_id")


## Biological processes

# Lets retreive the attributes for the marker set
GeneSet_annotated <- biomartr::biomart(genes = GEP44,
                                       mart       = "plants_mart",
                                       dataset    = "athaliana_eg_gene",
                                       attributes = attributes_to_retrieve,
                                       filters =  "ensembl_gene_id" )

# performing the over representation analysis (ORA) for Gene Ontology class - Biological processes
GeneSet_ORA_BP <- enrichGO(gene = GeneSet_annotated$entrezgene_id,
                           universe = species_genes_annotated$entrezgene_id,
                           OrgDb = org.At.tair.db,  # contains the TAIR/Ensembl id to GO correspondence for A. thaliana
                           keyType = "ENTREZID",
                           ont = "BP",
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05,
                           readable = TRUE,
                           pool = FALSE)
  

p <- dotplot(GeneSet_ORA_BP, showCategory = 10) + 
  labs(colour = "p.adj") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = grp_col[2]),
          axis.title.x = element_text(size = 36, face = "bold", color = "black"), 
          axis.ticks.length = unit(.30, "cm"), 
          axis.text.x = element_text(size = 36, color = "black", face = "bold"),
          axis.text.y = element_text(size = 36, color = "black", face = "bold"),
          legend.key = element_rect(size = 24),
          legend.text = element_text(size = 32),
          legend.title = element_text(size = 32, face = "bold"),
          legend.spacing = unit(2.0, 'cm'),
          legend.key.size = unit(3,"line"))
  
ggsave(filename = "Figure_1B.png", plot = p, width = 14, height = 14, dpi = 300)

