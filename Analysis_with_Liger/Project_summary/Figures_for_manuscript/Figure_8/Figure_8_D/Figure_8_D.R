# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the seurat object
with_l = loadRData("withCC.RData")

with_s = loadRData("seurat_withCC.RData")


cell_clusters = merge(with_l, with_s, by = "cell_ID")
rownames(cell_clusters) = cell_clusters$cell_ID

cell_clusters = cell_clusters[, c(2, 3)]

colnames(cell_clusters) = c("ID1", "ID2")

cell_clusters$ID1 <- factor(cell_clusters$ID1, levels = str_sort(levels(cell_clusters$ID1), numeric = TRUE))

cell_clusters$ID2 <- factor(cell_clusters$ID2, levels = str_sort(levels(cell_clusters$ID2), numeric = TRUE))


df <- cell_clusters %>%
  make_long(ID1, ID2)

df = df[order(as.numeric(df$node), decreasing = FALSE), ]

df$node = factor(df$node, levels = seq(0, 17, 1))

p <- ggplot(df, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = node, label = node)) +
  geom_sankey(flow.alpha = .7,
              node.color = "gray30") +
  geom_sankey_label(size = 10, color = "white", fill = "gray40") +
  theme_sankey(base_size = 18) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        axis.text = element_blank())

ggsave(plot = p, filename = "Trial2.png", width = 16, height = 12, dpi = 300)

p <- ggplot(df, aes(x = x, 
                    next_x = next_x, 
                    node = node, 
                    next_node = next_node,
                    fill = node, label = node)) +
  geom_sankey(flow.alpha = 1.2) +
  theme_sankey(base_size = 20) +
  labs(x = NULL) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),
        axis.text = element_blank())

ggsave(plot = p, filename = "Figure_8_D.png", width =14, height = 12, dpi = 300)

