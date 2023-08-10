gt <- ggplot_gtable(ggplot_build(p))

ge <- subset(gt$layout, name == "panel")

grid.draw(gt[ge$t:ge$b])