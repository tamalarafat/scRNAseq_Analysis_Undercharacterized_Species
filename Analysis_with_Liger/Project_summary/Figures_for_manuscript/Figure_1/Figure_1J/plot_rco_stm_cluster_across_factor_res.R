# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

###
# Load data tables
###

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/GET_RCO_STM/RCO_STM_Clustering.RData")

mdf$rco = 4
mdf$stm = 6

for(i in c(1:nrow(mdf))){
  if(mdf$stm_cluster[i] == mdf$rco_cluster[i]){
    mdf$rco[i] = 5
    mdf$stm[i] = 5
  }
}

mdf$K = as.factor(rownames(mdf))

mdf = mdf[, -c(1, 2)]
long_mat = mdf %>% melt(
  id.vars = "K",
  variable.name = "Cluster_type",
  value.name = "Clusters"
)

long_mat$Clusters = as.factor(long_mat$Clusters)

long_mat$Clusters = factor(long_mat$Clusters, levels = str_sort(levels(long_mat$Clusters), numeric = TRUE))



p <- ggplot(data = long_mat, aes(x = K, y = Clusters, group = Cluster_type)) + 
  geom_line(aes(color = Cluster_type), size = 2, position = position_dodge(width = 0.5)) + 
  scale_x_discrete(breaks = str_c("K_", seq(10, 50, 5)), labels = seq(10, 50, 5)) + 
  scale_color_manual(name = "Cluster", values = grp_col[c(3, 5)], labels = c("RCO", "STM")) + 
  scale_y_discrete(breaks = c(0, 5, 10, 15)) + 
  xlab("Factorization K") + 
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 48, colour = "black"), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 48, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_blank(),
    title = element_text(size = 40, face = "bold", colour = "black"),
    legend.text = element_text(size = 40, face = "bold", colour = "black"),
    legend.key.size = unit(1.5, 'cm'),
    legend.key.height= unit(0, 'cm'),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position = "top") + 
  guides(color = guide_legend(override.aes = list(size = 10)))

ggsave(filename = "Figure_1J.png", plot = p, width = 14, height = 14, dpi = 300)
