# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

###
# Load data tables
###


DF_OX_1E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/CH_01/Mean_Variance_CH_01.RData")
DF_OX_2E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/CH_02/Mean_Variance_CH_02.RData")
DF_OX_3E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/CH_03/Mean_Variance_CH_03.RData")
DF_OX_7E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/CH_07/Mean_Variance_CH_07.RData")

DF_COL_1E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/AT_01/Mean_Variance_AT_01.RData")
DF_COL_2E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/AT_02/Mean_Variance_AT_02.RData")
DF_COL_3E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/AT_03/Mean_Variance_AT_03.RData")
DF_COL_5E = loadRData("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/HVGs_for_each_dataset/Liger/3000/Ortho_datasets/AT_05/Mean_Variance_AT_05.RData")


a = cbind.data.frame(DF_OX_1E, DF_OX_2E, DF_OX_3E, DF_OX_7E, DF_COL_1E, DF_COL_2E, DF_COL_3E, DF_COL_5E)
grep(pattern = "avgExpression", colnames(a))

DFavg = rowMeans(a[, grep(pattern = "avgExpression", colnames(a))])
DFvar = rowMeans(a[, grep(pattern = "geneVariance", colnames(a))])

b = data.frame(avgExpression = DFavg, geneVariance = DFvar)
b <- b[!is.infinite(rowSums(b)),]

a = read_lines("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/Species_WT/Compare_HVGs_Liger/HVG_intersect_reps_union_species_without_mincells_2324.txt")
c = read_lines("/home/ytamal2/Documents/2022/Final_part_PhD_Thesis/Analyses/Species_WT/Compare_HVGs_Liger/Seurat_3000_HVG_intersect_reps_union_species_without_mincells_2294.txt")

# Define color
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A")

ggplot(b, aes(x = avgExpression, y = geneVariance)) + geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = b[a, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))

ggplot(b, aes(x = avgExpression, y = geneVariance)) + geom_point(alpha = 0.5, color = grp_col[2]) + 
  geom_point(data = b[c, ], aes(x = avgExpression, y = geneVariance), color = grp_col[3], alpha = 0.5) + 
  xlab("Gene Expression Mean (log10)") + 
  ylab("Gene Expression Variance (log10)") + 
  theme_classic() + 
  theme(
    line = element_blank(),
    panel.border = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 18, face = "bold"))
