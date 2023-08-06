# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/home/yasir/Documents/Chapter1/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

modules = str_sort(colnames(integrated.data@meta.data)[grep(x = colnames(integrated.data@meta.data), pattern = "Cluster")], numeric = TRUE)


for (j in c(1:length(modules))) {
  p <-
    FeaturePlot(
      integrated.data,
      features = modules[j],
      reduction = "umap",
      pt.size = 1.5,
      order = T,
      min.cutoff = 0.001,
      cols = c("grey", RColorBrewer::brewer.pal(9, "Reds")[9])
    ) +
    NoLegend() +
    theme(
      line = element_blank(),
      panel.border = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      plot.title = element_blank()
    )
  ggsave(
    plot = p,
    filename = str_c("Feature_plot_", modules[j], "_module_score.png"),
    width = 8,
    height = 8,
    dpi = 300
  )
}
