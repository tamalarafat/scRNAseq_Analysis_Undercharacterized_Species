# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Chapter1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
load("/home/yasir/Documents/Chapter1/UMAP_Selected/seurat_object_of_K_44_Without_GEP20_22.RData")

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

rownames(modulesAvg) <- str_c("cluster ", c(0:17))

colnames(modulesAvg) <- str_c("cluster ", c(0:17))


modulesAvg$clusters <- rownames(modulesAvg)

df_H <-
  modulesAvg %>% melt(
    id.vars = "clusters",
    variable.name = "modules",
    value.name = "score"
  )

ggplot(df_H, aes(clusters, modules)) +
  geom_tile(aes(fill = score), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")


cluster_cells = WhichCells(integrated.data, idents = "0")

a = integrated.data@meta.data[cluster_cells, modules]
b = data.frame(colMeans(a))
colnames(b)[1] = "cluster_0"
