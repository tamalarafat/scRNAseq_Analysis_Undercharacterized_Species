# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

# Save the annotation file containing Cardamine and Arabidopsis orthologues IDs with gene description
ATH = read.delim("/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Annotation_files/hiruta_version12_orthologs.txt", sep = "")

# Assign colnames of the columns of interest
colnames(ATH)[c(1, 2, 3, 4)] <- c("CH_ID", "AT_ID", "gene_length_CH", "gene_length_AT")

# Let's replace "_" with "-"
ATH$CH_ID = str_replace(ATH$CH_ID, pattern = "_", replacement = "-")

# Lets add RCO gene to the description file
rownames(ATH) <- ATH$CH_ID

# Adding the RCO orthologues Arabidopsis gene ID to the file
ATH["Chir06Ox-b35150.2", "AT_ID"] <- "AT5G67651"

rownames(ATH) <- NULL

write_csv(ATH, "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Beginning_of_a_compendium/Input_Data/Annotation_files/Gene_description_CH_V12_orthologs.csv")

