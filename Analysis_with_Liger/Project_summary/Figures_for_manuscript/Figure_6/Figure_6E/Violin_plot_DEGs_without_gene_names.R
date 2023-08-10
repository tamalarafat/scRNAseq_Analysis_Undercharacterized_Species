# Load all the functions stored in scripts from the folder housing the scripts
scripts_list <- list.files("/home/yasir/Documents/Thesis_PhD/Chapter_1/Functions", pattern = "*.R$", full.names = TRUE) 
sapply(scripts_list, source, .GlobalEnv)

## Define the colors
grp_col = c("#FD6467", "#00A08A", "#F98400", "#046C9A", "#075149FF", "#FFA319FF", "#00BF81", "#767676FF", "#FD8CC1FF")

# load the seurat object
load("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Seurat_object_without_GEP_20_22_CC/seurat_object_of_K_44_Without_GEP20_22_CC.RData")

Idents(integrated.data) <- "RNA_snn_res.0.3"

# Markers directory
markers_dir = "/home/yasir/Documents/Thesis_PhD/Chapter_1/Remove_GEP20_22_and_CCGEPs/Cluster_DEGs_Markers_Candidates/Cluster_biomarkers"

# Lets get the saved file names
temp_file_names = str_sort(list.files(path = markers_dir, pattern = "Cluster"), numeric = TRUE)

marker_list = list()

for (i in c(1:length(temp_file_names))){
  
  markers = read.csv(str_c(markers_dir, "/", temp_file_names[i]), row.names = 1)
  
  marker_list[[i]] <- rownames(markers)
  
  names(marker_list)[i] <- str_c("cluster", parse_number(temp_file_names[i]), "_module")
}

cluster_markers = data.frame(DEGs = c("Chir06Ox-b21610.2", "Chir02Ox-b11720.2", "Chir06Ox-b24210.2", "Chir01Ox-b41580.2", "Chir02Ox-b15590.2", "Chir08Ox-b23620.2", 
                                      "Chir06Ox-b35870.2", "Chir07Ox-b19660.2", "Chir03Ox-b28030.2", "Chir07Ox-b09460.2", "Chir01Ox-b32760.2", "Chir01Ox-b19700.2",
                                      "Chir05Ox-b05920.2", "Chir02Ox-b03100.2", "Chir03Ox-b35560.2", "Chir08Ox-b16630.2", "Chir06Ox-b13660.2", "Chir07Ox-b19200.2"))

# Load the gene description file
ATH = read.csv("/home/yasir/Documents/Thesis_PhD/Chapter_1/Analysis_objects/Annotation_files/Gene_description_CH_V12_orthologs.csv")
rownames(ATH) <- ATH$CH_ID

# Create ortho genes list and gene description file
description = ATH[rownames(ATH) %in% cluster_markers$DEGs, ]

# Description - All DEGs identified for the cluster
markers_description = merge(cluster_markers, description, by.x = "DEGs", by.y = "CH_ID", all.x = TRUE)
rownames(markers_description) = markers_description$DEGs

sorted_description = markers_description[cluster_markers$DEGs, ]

p1 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[1], pt.size = 0) + ylab(cluster_markers$DEGs[1]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p2 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[2], pt.size = 0) + ylab(cluster_markers$DEGs[2]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))


p3 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[3], pt.size = 0) + ylab(cluster_markers$DEGs[3]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p4 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[4], pt.size = 0) + ylab(cluster_markers$DEGs[4]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p5 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[5], pt.size = 0) + ylab(cluster_markers$DEGs[5]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p6 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[6], pt.size = 0) + ylab(cluster_markers$DEGs[6]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p7 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[7], pt.size = 0) + ylab(cluster_markers$DEGs[7]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p8 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[8], pt.size = 0) + ylab(cluster_markers$DEGs[8]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p9 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[9], pt.size = 0) + ylab(cluster_markers$DEGs[9]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p10 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[10], pt.size = 0) + ylab(cluster_markers$DEGs[10]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p11 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[11], pt.size = 0) + ylab(cluster_markers$DEGs[11]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p12 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[12], pt.size = 0) + ylab(cluster_markers$DEGs[12]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p13 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[13], pt.size = 0) + ylab(cluster_markers$DEGs[13]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p14 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[14], pt.size = 0) + ylab(cluster_markers$DEGs[14]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p15 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[15], pt.size = 0) + ylab(cluster_markers$DEGs[15]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p16 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[16], pt.size = 0) + ylab(cluster_markers$DEGs[16]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))

p17 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[17], pt.size = 0) + ylab(cluster_markers$DEGs[17]) + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "top",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))


p18 <- VlnPlot(integrated.data, features = cluster_markers$DEGs[18], pt.size = 0) + ylab(cluster_markers$DEGs[18]) + 
  xlab("Cell clusters") + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_blank(),
        legend.position = "none",
        legend.key.height = unit(0.8, 'cm'),
        legend.key.width = unit(0.8, 'cm'),
        legend.direction = "horizontal",
        legend.key = element_rect(size = 34),
        legend.spacing = unit(0.2, 'cm'),
        legend.text = element_text(size = 34, face = "bold"),
        legend.key.size = unit(4, "line")) + 
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 8)))


arranged_fig <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, ncol = 1, nrow = 18, legend = "none")
my_fig <- annotate_figure(arranged_fig, bottom = element_blank(), 
                          left = element_blank() + bgcolor("white"))
ggsave(
  "Violin_plot_DEGs_no_names_2.png",
  my_fig,
  width = 14,
  height = 20,
  dpi = 300
)

