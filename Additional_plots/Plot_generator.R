# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###############################

# Load the saved integrated file
integrated.data <- loadRData("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Analyses/Species_WT/Liger_analysis_30_09_2022/Liger_objects/Liger_object_K_50.RData")

integrated.data@clusters <- factor(integrated.data@clusters, levels = seq(1, length(levels(integrated.data@clusters))))

integrated.data@cell.data$dataset <- factor(integrated.data@cell.data$dataset, levels = c("WC1", "WC2", "WC3", "WC5", "WO1", "WO2", "WO3", "WO7"))

# What are the datasets
dfs <- levels(integrated.data@cell.data$dataset)

# Lets remove digits from strings, so we have only the species/genotypes initials
orgs <- unique(gsub("[0-9]+", "", dfs))

# if we have more than two organisms, we need to define that many number of empty variables
org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string
org2 = dfs[grep(pattern = orgs[2], dfs)] # output is a boolean for the whole string

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

thaliana_reps = df_list
hirsuta_reps = df_list
WO1_RS <- df_list[[1]]

p <- ggplot(WO1_RS, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  coord_cartesian(ylim = c(round(min(WO1_RS$MeanSpecificity)), round(max(WO1_RS$MeanSpecificity)))) +
  scale_y_continuous(breaks = seq(round(min(WO1_RS$MeanSpecificity)), round(max(WO1_RS$MeanSpecificity)), 10)) + 
  xlab("GEPs") + 
  ylab("Percent specificity (%)") + 
  coord_flip() + 
  theme(
    legend.key = element_rect(size = 20),
    axis.title.x = element_text(size = 38, face = "bold", colour = "black"), 
    axis.title.y = element_text(size = 38, face = "bold", colour = "black"), 
    axis.ticks.length = unit(.20, "cm"), 
    axis.text = element_text(size = 38, face = "bold", colour = "black"),
    title = element_text(size = 38, face = "bold", colour = "black"))

ggsave(filename = "Trial_WC1_replicate_specificity.png", plot = p, width = 24, height = 32, dpi = 300)
