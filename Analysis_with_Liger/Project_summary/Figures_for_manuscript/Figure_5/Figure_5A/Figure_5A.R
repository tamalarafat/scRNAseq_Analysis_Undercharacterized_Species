# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22/seurat_object_of_K_44_Without_GEP20_22.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Let's get the coefficient matrix
coef = as.data.frame(integrated.data@reductions[["inmfcc"]]@cell.embeddings)
colnames(coef) = str_c("GEP_", parse_number(colnames(coef)))

coef = order_factorized_matrix(coef)

cells_order = c()

for (i in levels(integrated.data)){
  cell_idents = WhichCells(integrated.data, idents = i)
  cells_order = c(cells_order, cell_idents)
}

mat_H = coef[cells_order, ]

mat_H = data.frame(Cells = rownames(mat_H), mat_H)
mat_H$Cells = factor(mat_H$Cells, levels = mat_H$Cells)


mat_H <-
  mat_H %>% melt(
    id.vars = "Cells",
    variable.name = "GEPs",
    value.name = "Expression"
  )

intercepts = c()

for (i in levels(integrated.data)) {intercepts = c(intercepts, length(WhichCells(integrated.data, idents = i)))}

intercepts = cumsum(intercepts)[-length(intercepts)]

p <-
  ggplot(mat_H, aes(x = Cells, y = GEPs)) + geom_tile(aes(fill = Expression)) +
  scale_fill_gradient(
    name = "Rest\nof the\ncells",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  new_scale("fill") +
  geom_tile(aes(fill = Expression), data = mat_H[which(mat_H$GEPs %in% str_c("GEP_", c(8, 11, 25, 44))),]) +
  scale_fill_gradient(
    name = "cells with\ncell-cyle\nactivity",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Reds")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  geom_vline(xintercept = intercepts, linetype = "dotted", size = 1) + 
  xlab("Cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(21, seq(2, 44, 2))), labels = c(21, seq(2, 44, 2))) + 
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
ggsave(
  filename = "Figure_5A_GEPs_marked.png",
  plot = p,
  width = 28,
  height = 18,
  dpi = 300
)

# No legend
p <-
  ggplot(mat_H, aes(x = Cells, y = GEPs)) + geom_tile(aes(fill = Expression)) +
  scale_fill_gradient(
    name = "Rest\nof the\ncells",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Blues")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  new_scale("fill") +
  geom_tile(aes(fill = Expression), data = mat_H[which(mat_H$GEPs %in% str_c("GEP_", c(8, 11, 25, 44))),]) +
  scale_fill_gradient(
    name = "cells with\ncell-cyle\nactivity",
    low = "white",
    high = RColorBrewer::brewer.pal(9, "Reds")[8],
    limit = c(min(mat_H$Expression), 1),
    space = "Lab",
    guide = "colourbar"
  ) +
  geom_vline(xintercept = intercepts, linetype = "dotted", size = 1) + 
  xlab("Cells") +
  scale_y_discrete(breaks = str_c("GEP_", c(21, seq(2, 44, 2))), labels = c(21, seq(2, 44, 2))) + 
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
ggsave(
  filename = "Figure_5A_GEPs_marked_2.png",
  plot = p,
  width = 28,
  height = 18,
  dpi = 300
)
