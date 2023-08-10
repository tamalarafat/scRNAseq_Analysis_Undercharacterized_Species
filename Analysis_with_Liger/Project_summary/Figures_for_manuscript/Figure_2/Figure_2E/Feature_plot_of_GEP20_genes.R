# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF")

load("/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/On_coefficient/Seurat_objects/seurat_object_of_K_44.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"


GEP_genes = read.delim("GEP_20_genes.txt", header = FALSE, col.names = "gene_ID")
GEP_genes = GEP_genes$gene_ID



for (j in c(1:length(GEP_genes))){
  
  base <- FeaturePlot(integrated.data, features = GEP_genes[j], reduction = "umap", pt.size = 1, min.cutoff = 0.001, split.by = "Replicates", order = T, combine = FALSE, cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9]))
  
  temp_list = list()
  
  for (i in c(1:length(base))){
    p <- base[[i]]
    p <- p + NoLegend() + 
      labs(title = str_replace(string = names(which(table(p[["data"]][["ident"]]) != 0)), pattern = "-", replacement = " ")) + 
      theme(
        axis.line = element_blank(),
        axis.title.y.right = element_blank(),
        panel.border = element_blank(),
        line = element_blank(), 
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(size = 28, face = "italic")) +  
      scale_color_gradient(low = "gray", high = RColorBrewer::brewer.pal(9, "Reds")[9], limits = c(0, max(integrated.data@assays$RNA@data[GEP_genes[j], ])), name = GEP_genes[j])
  
    temp_list[[i]] <- p
    }
  
  final_fig <- ggarrange(plotlist = temp_list, ncol = 4, nrow = 1)
  
  ggsave(filename = str_c(GEP_genes[j], "_expression_across_reps.png"), plot = final_fig, 
         width = 14, height = 5, dpi = 300, bg = "white")
}


