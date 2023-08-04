# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Calculate alignment and agreement score for all factorizaation rank K = 10 - 50
###

# Directory containing all the Liger objects
Liger_objects_dir = "/biodata/dep_tsiantis/grp_laurent/tamal/2023/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/Liger_objects"

# Directory to store the output
storing_dir = "/netscratch/dep_tsiantis/grp_laurent/tamal/2023/Analyses/Analysis_of_single_species_Cardamine/Liger_Analyses/Liger_analysis_strategy_1/Alignment_and_agreement_scores"

# Applying the function to calculate alignment and agreement for each factorization rank K
calc_alignment_and_agreement(objects_dir = Liger_objects_dir, 
                             file_name_pattern = "Liger_object_K", 
                             store_dir = storing_dir)
