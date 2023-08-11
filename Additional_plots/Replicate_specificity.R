# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###############################

##### Functions used in this script

DS_Groups <- function(LigerObject){
  # What are the datasets
  dfs <- levels(LigerObject@cell.data$dataset)
  
  # Lets remove digits from strings, so we have only the species/genotypes initials
  orgs <- unique(gsub("[0-9]+", "", dfs))
  
  # if we have more than two organisms, we need to define that many number of empty variables
  org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string
  org2 = dfs[grep(pattern = orgs[2], dfs)] # output is a boolean for the whole string
  
  # Let's create a directory to store the figures
  
  for (i in c(1:length(org1))){
    for (j in c(1:length(org2))){
      DS <- calcDatasetSpecificity(LigerObject, dataset1 = org1[i], dataset2 = org2[j], do.plot = F)
      
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
  }
  DF$MeanSpecificity = rowMeans(DF[, -c(1)])
  DF = DF[ ,c("Factors", "MeanSpecificity")]
  return(DF)
} 

##### Here loading all the packages and functions for this script to run ######
# Theme specification
themeSpec <- theme(axis.title.x = element_text(size = 34, face = "bold"), 
                   axis.title.y = element_text(size = 34, face = "bold"), 
                   axis.ticks.length = unit(.30, "cm"), 
                   axis.text = element_text(size = 34, face = "bold"),
                   legend.key.size = unit(4,"line"),
                   legend.key = element_rect(size = 34),
                   legend.text = element_text(size = 26, face = "bold"))

theme_set(theme_bw())

###############################

# Load the saved integrated file
integrated.data <- loadRData("/netscratch/dep_tsiantis/grp_laurent/tamal/2022/Analyses/Species_WT/Liger_analysis_30_09_2022/Liger_objects/Liger_object_K_50.RData")

integrated.data@clusters <- factor(integrated.data@clusters, levels = seq(1, length(levels(integrated.data@clusters))))

integrated.data@cell.data$dataset <- factor(integrated.data@cell.data$dataset, levels = c("WC1", "WC2", "WC3", "WC5", "WO1", "WO2", "WO3", "WO7"))

# Which genes to plot
HVGs <- integrated.data@var.genes

####### Generate mean species-specificity score

df_specificity = DS_Groups(integrated.data)

p <- ggplot(df_specificity, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  coord_cartesian(ylim = c(round(min(df_specificity$MeanSpecificity)), round(max(df_specificity$MeanSpecificity)))) +
  scale_y_continuous(breaks = seq(round(min(df_specificity$MeanSpecificity)), round(max(df_specificity$MeanSpecificity)), 10)) + 
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

ggsave(filename = "Average_replicate_specificity.png", plot = p, width = 24, height = 32, dpi = 300)

