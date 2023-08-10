# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the gene cluster file that is stored in the basis matrix dataframe
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Basis_objects/Basis.RData")

Basis_subset = Basis[, str_c("K_", c(10, 20, 30, 44, 50))]

p <- clustree(Basis_subset, prefix = "K_", node_text_size = 0, prop_filter = 0.2) + 
  labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
  theme(legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 34, face = "bold", color = "black"),
        legend.title = element_text(size = 34, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size = 12)))

ggsave(filename = "Figure_1A.png", plot = p, width = 13, height = 13, dpi = 600)



p <- clustree(Basis_subset, prefix = "K_", node_text_size = 10, prop_filter = 0.2) + 
  labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
  theme(legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 34, face = "bold", color = "black"),
        legend.title = element_text(size = 34, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size = 12)))

ggsave(filename = "Figure_1A_2.png", plot = p, width = 14, height = 13, dpi = 600)


# no legend
p <- clustree(Basis_subset, prefix = "K_", node_text_size = 0, prop_filter = 0.2) + 
  labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + NoLegend()
  

ggsave(filename = "Figure_1A_3.png", plot = p, width = 10, height = 10, dpi = 600)

# legend
p <- clustree(Basis_subset, prefix = "K_", node_text_size = 0, prop_filter = 0.2) + 
  labs(colour = "Factorization\nK", edge_colour = "Count", edge_alpha = "Proportion") + 
  theme(legend.key.size = unit(2,"line"),
        legend.text = element_text(size = 48, color = "black"),
        legend.title = element_text(size = 48, face = "bold", color = "black")) + 
  guides(colour = guide_legend(override.aes = list(size = 12)))

ggsave(filename = "Figure_1A_4.png", plot = p, width = 20, height = 20, dpi = 300)

