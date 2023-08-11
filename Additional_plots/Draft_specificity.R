# What are the datasets
dfs <- levels(integrated.data@cell.data$dataset)

# Lets remove digits from strings, so we have only the species/genotypes initials
orgs <- unique(gsub("[0-9]+", "", dfs))

# if we have more than two organisms, we need to define that many number of empty variables
org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string
org2 = dfs[grep(pattern = orgs[2], dfs)] # output is a boolean for the whole string

ml = list()

# Let's create a directory to store the figures

for (j in c(1:length(org2))){
    DS <- calcDatasetSpecificity(integrated.data, dataset1 = org1[1], dataset2 = org2[j], do.plot = F)
    
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

ml[[org1[1]]] <- DF

DS <- calcDatasetSpecificity(integrated.data, dataset1 = org1[1], dataset2 = org1[2], do.plot = F)
DF <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
DF <- DF[ ,c(2,1)]
DF$Factors <- factor(DF$Factors, levels = seq(1, length(DF$Factors)))

DS2 <- calcDatasetSpecificity(integrated.data, dataset1 = org1[1], dataset2 = org1[4], do.plot = F)
temp_DF <- data.frame(Specificity = DS2[[3]], Factors = c(1:length(DS2[[3]])))
temp_DF <- temp_DF[ ,c(2,1)]
temp_DF$Factors <- factor(temp_DF$Factors, levels = seq(1, length(temp_DF$Factors)))

DF <- cbind.data.frame(DF, temp_DF$Specificity)

DF$MeanSpecificity = rowMeans(DF[, -c(1)])
DF = DF[ ,c("Factors", "MeanSpecificity")]

p <- ggplot(DF, aes(x = Factors, y = MeanSpecificity)) + 
  geom_bar(fill = "#075149FF", stat = "identity") + 
  coord_cartesian(ylim = c(round(min(DF$MeanSpecificity)), round(max(DF$MeanSpecificity)))) +
  scale_y_continuous(breaks = seq(round(min(DF$MeanSpecificity)), round(max(DF$MeanSpecificity)), 10)) + 
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

ggsave(filename = "Trial_WC1.png", plot = p, width = 24, height = 32, dpi = 300)

DS <- calcDatasetSpecificity(integrated.data, dataset1 = org1[1], dataset2 = org1[2], do.plot = F)
DF <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
DF <- DF[ ,c(2,1)]
DF$Factors <- factor(DF$Factors, levels = seq(1, length(DF$Factors)))

DS <- calcDatasetSpecificity(integrated.data, dataset1 = org1[1], dataset2 = org1[2], do.plot = F)
DF <- data.frame(Specificity = DS[[3]], Factors = c(1:length(DS[[3]])))
DF <- DF[ ,c(2,1)]
DF$Factors <- factor(DF$Factors, levels = seq(1, length(DF$Factors)))
# Here let's create one replicate vs all replicates of a species
# So the comparison will be meana specificity of a replicate to species specificity

DS_replicates<- function(LigerObject){
  # What are the datasets
  dfs <- levels(LigerObject@cell.data$dataset)
  
  # Lets remove digits from strings, so we have only the species/genotypes initials
  orgs <- unique(gsub("[0-9]+", "", dfs))
  
  # if we have more than two organisms, we need to define that many number of empty variables
  org1 = dfs[str_detect(dfs, pattern = orgs[1])] # output is a boolean for the whole string
  org2 = dfs[grep(pattern = orgs[2], dfs)] # output is a boolean for the whole string
  
  # Let's create an empty list to store the specificity scores
  df_list = list()
  
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
  ml[[org1[1]]] <- DF
  }
  DF$MeanSpecificity = rowMeans(DF[, -c(1)])
  DF = DF[ ,c("Factors", "MeanSpecificity")]
  return(DF)
} 