# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

modules = str_sort(colnames(integrated.data@meta.data)[grep(x = colnames(integrated.data@meta.data), pattern = "Cluster")], numeric = TRUE)

clusterID = str_sort(levels(Idents(integrated.data)), numeric = TRUE)

avg_list = list()

for(i in c(1:length(clusterID))){
  cluster_cells = WhichCells(integrated.data, idents = clusterID[i])
  
  df = integrated.data@meta.data[cluster_cells, modules]
  
  avgdf = data.frame(colMeans(df))
  
  colnames(avgdf) <- str_c("cluster_", clusterID[i])
  
  avg_list[[i]] = avgdf
}

modulesAvg = do.call(cbind.data.frame, avg_list)

colnames(modulesAvg) <- str_c("biomarkers_cluster", c(0:17))

rownames(modulesAvg) <- str_c("cluster", c(0:17))

modulesAvg$clusters <- rownames(modulesAvg)

df_H <-
  modulesAvg %>% melt(
    id.vars = "clusters",
    variable.name = "modules",
    value.name = "score"
  )

df_H$clusters = factor(df_H$clusters, levels = str_c("cluster", c(0:17)))

df_H$modules = factor(df_H$modules, levels = str_c("biomarkers_cluster", c(0:17)))

p <- ggplot(df_H, aes(modules, clusters)) +
  geom_tile(aes(fill = score), colour = "white") +
  scale_fill_gradient(
    name = "Score",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$score), max(df_H$score)),
    space = "Lab",
    guide = "colourbar"
  ) + 
  scale_x_discrete(breaks = str_c("biomarkers_cluster", c(0:17)), labels = str_c("set ", c(0:17))) + 
  scale_y_discrete(breaks = str_c("cluster", c(0:17)), labels = c(0:17)) + 
  xlab("Cell-type markers") + 
  ylab("Cell clusters") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 48, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 48, colour = "black", vjust = 0.5),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  )
ggsave(plot = p, filename = "Figure_6_C.png", width = 14, height = 14, dpi = 300)

# No legend
p <- ggplot(df_H, aes(modules, clusters)) +
  geom_tile(aes(fill = score), colour = "white") +
  scale_fill_gradient(
    name = "Score",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(df_H$score), max(df_H$score)),
    space = "Lab",
    guide = "colourbar"
  ) + 
  scale_x_discrete(breaks = str_c("biomarkers_cluster", c(0:17)), labels = str_c("set ", c(0:17))) + 
  scale_y_discrete(breaks = str_c("cluster", c(0:17)), labels = c(0:17)) + 
  xlab("Cell-type markers") + 
  ylab("Cell clusters") + 
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 48, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 48, colour = "black", vjust = 0.5),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  ) + NoLegend()
ggsave(plot = p, filename = "Figure_6_C2.png", width = 14, height = 14, dpi = 300)

