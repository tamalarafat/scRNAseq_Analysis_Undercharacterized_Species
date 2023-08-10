# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

# Get the coefficient matrix
coef = as.data.frame(integrated.data@reductions[["inmf"]]@cell.embeddings)
colnames(coef) = str_c("GEP_", c(1:44))

coef = order_factorized_matrix(coef)

mat_H = data.frame(Cells = rownames(coef), coef)
mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)

mat_H <-
  mat_H %>% melt(
    id.vars = "Cells",
    variable.name = "GEPs",
    value.name = "Expression"
  )

p <-
  ggplot(mat_H, aes(x = Cells, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion\nusage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  scale_y_discrete(breaks = str_c("GEP_", seq(2, 44, 2)), labels = seq(2, 44, 2)) + 
  xlab("Cells") +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, colour = "black"),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  )

ggsave(filename = "Figure_1F.png", plot = p, width = 28, height = 18, dpi = 300)


# No legend
p <-
  ggplot(mat_H, aes(x = Cells, y = GEPs, fill = Expression)) + geom_tile() +
  scale_fill_gradient(
    name = "Proportion\nusage",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  scale_y_discrete(breaks = str_c("GEP_", seq(2, 44, 2)), labels = seq(2, 44, 2)) + 
  xlab("Cells") +
  theme(
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.title = element_text(size = 48, color = "black"),
    axis.ticks.length = unit(.20, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 48, colour = "black"),
    legend.title = element_text(size = 42, colour = "black"),
    legend.key.size = unit(3, "line"),
    legend.text = element_text(size = 48)
  ) + NoLegend()

ggsave(filename = "Figure_1F_2.png", plot = p, width = 28, height = 18, dpi = 300)

