# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object/seurat_object_of_K_44.RData")

###
# Load data tables
###

load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Suggest_K_and_Lambda/Suggest_K/suggest_K_ox_wt_liger.Rdata")

df = sg_K$data
df = df[df$calc == "KL_div", ]


p <- ggplot(df, aes(x = k, y = median_kl)) + 
  geom_line(color = grp_col[3], size = 6) + 
  geom_point(color = grp_col[4], size = 6) + 
  xlab("K (GEPs)") + 
  ylab("Median KL divergence") + 
  scale_x_continuous(breaks = seq(10, 60, 10)) + 
  theme(
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2], size = 2),
    axis.title = element_text(size = 48, colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 48, colour = "black"),
    title = element_text(size = 48, colour = "black"))

ggsave(filename = "Figure_1H.png", plot = p, width = 12, height = 12, dpi = 300)

# Sources
# https://github.com/welch-lab/liger/issues/103
# https://github.com/welch-lab/liger/issues/116
# https://www.rdocumentation.org/packages/rliger/versions/0.5.0/topics/suggestK
# https://broadinstitute.github.io/2020_scWorkshop/batch-correction-lab.html

