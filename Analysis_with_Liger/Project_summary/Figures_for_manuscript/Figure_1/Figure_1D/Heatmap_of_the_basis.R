# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Load the Seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

basis = as.data.frame(integrated.data@reductions[["inmf"]]@feature.loadings)
colnames(basis) = str_c("GEP_", c(1:44))

basis = order_factorized_matrix(basis)

mat_W = data.frame(Genes = rownames(basis), basis)
mat_W$Genes <- factor(mat_W$Genes, levels = mat_W$Genes)
mat_W <-
  mat_W %>% melt(
    id.vars = "Genes",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <- ggplot(mat_W, aes(x = GEPs, y = Genes, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Coefficient\nvalue",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_W$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("GEPs") +
  scale_x_discrete(breaks = str_c("GEP_", seq(2, 44, 2)), labels = seq(2, 44, 2)) + 
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
    axis.text.y = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  )
ggsave(filename = "Figure_1D.png", plot = p, width = 16, height = 18, dpi = 300)


# Without legend
p <- ggplot(mat_W, aes(x = GEPs, y = Genes, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Coefficient\nvalue",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_W$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  xlab("GEPs") +
  scale_x_discrete(breaks = str_c("GEP_", seq(2, 44, 2)), labels = seq(2, 44, 2)) + 
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
    axis.text.y = element_blank(),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  ) + NoLegend()
ggsave(filename = "Figure_1D_2.png", plot = p, width = 16, height = 18, dpi = 300)
