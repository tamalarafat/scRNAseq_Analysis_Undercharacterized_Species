# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# Replicate specificity
###############################

# Load the saved integrated file
integrated.data <- loadRData("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Liger_object/Liger_object_K_44.RData")

# Which genes to plot
HVGs <- integrated.data@var.genes

integrated.data@clusters <- factor(integrated.data@clusters, levels = seq(1, length(levels(integrated.data@clusters))))

integrated.data@cell.data$dataset <- factor(integrated.data@cell.data$dataset, levels = c("WO1", "WO2", "WO3", "WO7"))

# What are the datasets
dfs <- levels(integrated.data@cell.data$dataset)

# Lets remove digits from strings, so we have only the species/genotypes initials
orgs <- unique(gsub("[0-9]+", "", dfs))

# if we have more than two organisms, we need to define that many number of empty variables
org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string

# Let's create a directory to store the figures

# Let's create an empty list to store the specificity scores
df_list = list()

for (i in c(1:length(org1))){
  target_replicate = org1[i]
  remaining_replicates = org1[-c(i)]
  
  for (j in c(1:length(remaining_replicates))){
    DS <- calcDatasetSpecificity(integrated.data, dataset1 = target_replicate, dataset2 = remaining_replicates[j], do.plot = F)
    
    if (!exists("DF")){
      DF <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
      DF <- DF[ ,c(2,1)]
      DF$Factors <- factor(DF$Factors, levels = seq(1, length(DF$Factors)))
    }
    else {
      temp_df <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
      DF <- cbind.data.frame(DF, temp_df$Specificity)
    }
  }
  DF$MeanSpecificity = rowMeans(DF[, -c(1)])
  DF = DF[ ,c("Factors", "MeanSpecificity")]
  df_list[[org1[i]]] <- DF
  rm(DF)
}

hirsuta_reps = df_list

# OX 1
WO1_RS <- df_list[[1]]

p <- ggplot(WO1_RS, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 10), limits = c(-20, 20)) + 
  scale_x_discrete(breaks = seq(2, 44, 2), labels = seq(2, 44, 2)) + 
  coord_flip() +
  xlab("GEPs") + 
  ylab("Percent specificity (%)") + 
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 40, face = "bold", colour = "black"), 
    axis.title.y = element_text(size = 40, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 40, colour = "black"),
    title = element_text(size = 38, face = "bold", colour = "black"))

ggsave(filename = "Figure_1E_OX1.png", plot = p, width = 10, height = 14, dpi = 300, bg = "white")

# OX 2
WO2_RS <- df_list[[2]]

p1 <- ggplot(WO2_RS, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 10), limits = c(-20, 20)) + 
  scale_x_discrete(breaks = seq(2, 44, 2), labels = seq(2, 44, 2)) + 
  coord_flip() +
  xlab("GEPs") + 
  ylab("Percent specificity (%)") + 
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 40, face = "bold", colour = "black"), 
    axis.title.y = element_text(size = 40, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 40, colour = "black"),
    title = element_text(size = 38, face = "bold", colour = "black"))

ggsave(filename = "Figure_1E_OX2.png", plot = p1, width = 10, height = 14, dpi = 300, bg = "white")

# OX 3
WO3_RS <- df_list[[3]]

p2 <- ggplot(WO3_RS, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 10), limits = c(-20, 20)) + 
  scale_x_discrete(breaks = seq(2, 44, 2), labels = seq(2, 44, 2)) + 
  coord_flip() +
  xlab("GEPs") + 
  ylab("Percent specificity (%)") + 
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 40, face = "bold", colour = "black"), 
    axis.title.y = element_text(size = 40, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 40, colour = "black"),
    title = element_text(size = 38, face = "bold", colour = "black"))

ggsave(filename = "Figure_1E_OX3.png", plot = p2, width = 10, height = 14, dpi = 300, bg = "white")

# OX 7
WO7_RS <- df_list[[4]]

p3 <- ggplot(WO7_RS, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  scale_y_continuous(breaks = seq(from = -20, to = 20, by = 10), limits = c(-20, 20)) + 
  scale_x_discrete(breaks = seq(2, 44, 2), labels = seq(2, 44, 2)) + 
  coord_flip() +
  xlab("GEPs") + 
  ylab("Percent specificity (%)") + 
  theme(
    panel.border = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = grp_col[2]),
    axis.title.x = element_text(size = 40, face = "bold", colour = "black"), 
    axis.title.y = element_text(size = 40, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 40, colour = "black"),
    title = element_text(size = 38, face = "bold", colour = "black"))

ggsave(filename = "Figure_1E_OX7.png", plot = p3, width = 10, height = 14, dpi = 300, bg = "white")

