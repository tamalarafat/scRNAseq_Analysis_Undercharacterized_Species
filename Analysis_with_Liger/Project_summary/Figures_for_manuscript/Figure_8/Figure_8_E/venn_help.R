# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(x, fill = c(grp_col[2], grp_col[1], grp_col[4], grp_col[9]), # Circles
             lwd = 2,
             lty = 'blank',
             # Numbers
             cex = 5,
             # Set names
             cat.cex = 1,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.dist = c(0.055, 0.055, 0.1, 0.1))