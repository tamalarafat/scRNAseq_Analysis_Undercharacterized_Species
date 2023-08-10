# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

df  = loadRData("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Alignment_and_agreement_scores/alignment_and_agreement_score.RData")

df$GEP = parse_number(rownames(df))
# df = df[df$GEP %in% seq(10, 50, 5), ]

p <- ggplot(df, aes(x = GEP)) + 
  geom_line(aes(y = Alignment, colour = "Alignment"), size = 6) + 
  geom_line(aes(y = Agreement, colour = "Agreement"), size = 6) + 
  geom_point(aes(y = Alignment), color = grp_col[4], size = 2) + 
  geom_point(aes(y = Agreement), color = grp_col[4], size = 2) + 
  xlab("K (GEPs)") + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual("Score", breaks = c("Alignment", "Agreement"),
                      values = c("Alignment" = grp_col[3], "Agreement" = grp_col[1])) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 48, colour = "black"), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 48, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 48, colour = "black"),
    title = element_text(size = 48, colour = "black"),
    legend.text = element_text(size = 48, colour = "black"),
    legend.key.size = unit(1, 'cm'),
    legend.key.height= unit(2, 'cm'),
    legend.key = element_rect(colour = "transparent", fill = "white"),
    legend.position = "top")

ggsave(filename = "Figure_1C.png", plot = p, width = 12, height = 12, dpi = 300)


p <- ggplot(df, aes(x = GEP)) + 
  geom_line(aes(y = Alignment, colour = "Alignment"), size = 6) + 
  geom_line(aes(y = Agreement, colour = "Agreement"), size = 6) + 
  geom_point(aes(y = Alignment), color = grp_col[4], size = 2) + 
  geom_point(aes(y = Agreement), color = grp_col[4], size = 2) + 
  xlab("K (GEPs)") + 
  scale_x_continuous(breaks = seq(10, 50, 5), labels = str_c("GEP", seq(10, 50, 5))) + 
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual("Score", breaks = c("Alignment", "Agreement"),
                      values = c("Alignment" = grp_col[3], "Agreement" = grp_col[1])) + 
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 48, colour = "black"), 
    axis.title.y = element_blank(), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text.x = element_text(size = 48, colour = "black", angle = 90, vjust = 0.5),
    axis.text.y = element_text(size = 48, colour = "black"),
    title = element_text(size = 48, colour = "black"),
    legend.text = element_text(size = 48, colour = "black"),
    legend.key.size = unit(1, 'cm'),
    legend.key.height= unit(2, 'cm'),
    legend.key = element_rect(colour = "transparent", fill = "white"))

ggsave(filename = "Figure_1C_2.png", plot = p, width = 12, height = 12, dpi = 300)

