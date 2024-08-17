library(dplyr)
library(Seurat)
library(MASS)
library(ggplot2)
library(caret)
library(Metrics)
library(viridis)
library(ggbeeswarm)

load("WT.Robj")

# Add annotations to WT object
cluster_order <- c(2, 38, 6, 53, 41, 10, 48, 17, 32, 37, 13, 52, 3, 44, 21, 34, 14, 45, 46, 50, 28, 4, 42, 54, 0, 30, 16, 47, 20, 18, 25, 33, 7, 9, 24, 31, 8, 49, 51, 43, 29, 5, 19, 23, 35, 26, 22, 39, 27, 11, 15, 1, 40, 36, 12)
cluster_names <- c('2_Epiblast', '38_PostProxPrimStreak', '6_EarlyMeso', '53_PGC', '41_BodyWall', '10_ProxExE', '48_DiffExE', '17_DiffExE', '32_Trophoblasts', '37_DefEndo', '13_PrimEndo', '52_YolkSac2', '3_YolkSac1', '44_ParietalEndo', '21_GutEndo', '34_HemogenEndothel', '14_Angioblast', '45_MonocMacrophProg', '46_EMP', '50_HematoEndoProg', '28_PrimBloodProg', '4_Blood_early', '42_PrimFetal_Blood', '54_PrimFetal_Blood', '0_Blood_late', '30_Xmeso_early', '16_Allantois', '47_AmnioticMeso3', '20_AmnioticMeso2', '18_AmnioticMeso1', '25_PresomMeso', '33_EarlyMeso_post', '7_PostLat_IntermMeso', '9_SecHF', '24_PharyArchMeso', '31_PrimHF', '8_Somites', '49_GenitourMeso', '51_Artifact', '43_NodeNotochord', '29_IndEpi_early', '5_IndEpi_late', '19_SurfaceEcto', '23_NeuralRidge_post', '35_NeuralRidge_ant', '26_PrimStreak', '22_NMP', '39_PreneuralPlate_post', '27_PreneuralPlate_ant', '11_Forebrain', '15_Midbrain', '1_NeuralPlate', '40_FloorPlate', '36_MotorNeuron', '12_NeuralCrest')
cluster_colors <- c('#000000', '#474747', '#a0a0a0', '#DB083D', '#722c61', '#fcb7b7', '#f4595c', '#db484a', '#5d0000', '#f6c100', '#f7a800', '#f89e00', '#f98e00', '#fb7100', '#f5f200', '#fad4a6', '#deb88c', '#cba47a', '#bf986f', '#b28b63', '#a57d57', '#98704b', '#8b633f', '#7c5432', '#673f1f', '#2cfeff', '#5bdffe', '#7acafe', '#a9abfd', '#da8afc', '#d767ff', '#a746d0', '#9339bd', '#7927a4', '#62188e', '#3f006c', '#166788', '#099be2', '#492b00', '#1142fc', '#c2ecb3', '#99e2a9', '#4ad094', '#25c88a', '#00bf80', '#00ffda', '#00efca', '#00e1bc', '#00d2ac', '#00c19a', '#00ac83', '#009e75', '#008d62', '#00784b', '#003d0e')
cluster_group <- c('Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Placenta', 'Placenta', 'Placenta', 'Placenta', 'Embryo', 'YolkSac', 'YolkSac', 'YolkSac', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo')
cluster_germlayer <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'Meso', 'XEcto', 'XEcto', 'XEcto', 'XEcto', 'Endo', 'XEndo', 'XEndo', 'XEndo', 'XEndo', 'Endo', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'XEndo', 'Meso', 'Epiblast', 'Epiblast', 'Ecto', 'Ecto', 'Ecto', 'Meso', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto')

annotation <- data.frame(State=cluster_order, cluster_colors, cluster_names, cluster_group, cluster_germlayer)
rownames(annotation) <- annotation$cluster_names
annotation$State <- as.character(annotation$State)
merged_data <- WT@meta.data %>% left_join(annotation, by= c("State"="State"))

WT <- AddMetaData(WT, metadata = merged_data$cluster_germlayer, col.name = "cluster_germlayer")
WT <- AddMetaData(WT, metadata = merged_data$cluster_group, col.name = "cluster_group")
WT <- AddMetaData(WT, metadata = merged_data$cluster_names, col.name = "cluster_names")
WT <- AddMetaData(WT, metadata = merged_data$cluster_colors, col.name = "cluster_colors")

# All embryonic stages 
stages <- unique(WT@meta.data$stage)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, stages)
named_colors <- setNames(cluster_colors, cluster_names)
stages_numeric <- c(6.5, 7.0, 7.5, 8.0, 8.5, 9.0)

# Define a set of distinct colors for each stage (same as in PCA)
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames(stage_colors, e_days)
color_palette <- setNames(stage_colors, c("6.5","7","7.5","8","8.5","9"))
color_mapping2 <- setNames( stage_colors, c("6.5","7.0","7.5","8.0","8.5","9.0"))

# Apply LDA on Knockouts:
load("Eed_E85.Robj")
load("Rnf2_E85.Robj")
load("Kdm6a_E85.Robj")
load("Kdm6ayb_E85.Robj")

ko_objs <- c("Eed85","Rnf285", "Kdm6a85", "Kdm6ayb85")

for(ko_obj_name in ko_objs){
  ko_obj <- get(ko_obj_name)
  merged_data <- ko_obj@meta.data %>% left_join(annotation, by= c("State"="State"))

  ko_obj <- AddMetaData(ko_obj, metadata = merged_data$cluster_germlayer, col.name = "cluster_germlayer")
  ko_obj <- AddMetaData(ko_obj, metadata = merged_data$cluster_group, col.name = "cluster_group")
  ko_obj <- AddMetaData(ko_obj, metadata = merged_data$cluster_names, col.name = "cluster_names")
  ko_obj <- AddMetaData(ko_obj, metadata = merged_data$cluster_colors, col.name = "cluster_colors")
  assign(ko_obj_name, ko_obj)
}
###########################################################################################################
# use whole WT for training LDA: 
###########################################################################################################
train <- WT
Idents(train) <- train$stage
all_markers <- list()

for(s in stages) {
  # Find marker genes for cells of a specific stage compared to all others (only positive ones)
  markers <- FindMarkers(train, only.pos=TRUE, ident.1 = s, min.pct = 0.05, logfc.threshold = 0.25)
  all_markers[[s]] <- markers
}
lapply(all_markers, head)

# Filter markers based on p-value threshold (p < 0.05) for each stage
all_markers <- lapply(all_markers, function(markers) {
  markers[markers$p_val < 0.05, ]
})

# Create a vector of unique markers from all stages
selected_genes <- unique(unlist(lapply(all_markers, rownames)))
length(selected_genes) #should be 877
write.csv(data.frame(Gene = selected_genes), "pos_marker_genes_WT.csv", row.names = FALSE)
#selected_genes <- read.csv("pos_marker_genes_WT.csv", header=TRUE)$Gene

# Subset to include only selected genes
train_subset <- subset(train, features = selected_genes)

# Extract normalized expression data 
data_for_lda <- GetAssayData(train_subset, slot = "data")
data_for_lda <- as.matrix(data_for_lda)
data_for_lda <- t(data_for_lda)  

class_labels <- train_subset@meta.data$stage

# Running LDA on WT using uniform distributed prior
lda_result <- lda(data_for_lda, prior=rep(1/6,6), grouping = class_labels)
save(lda_result, file = "lda_result_WT.RData")
#load("lda_result_WT.RData")

lda_prediction <- predict(lda_result, data_for_lda)
str(lda_prediction)
lda_scales <- as.data.frame(lda_prediction$x)
lda_scales$stage <- class_labels

# abs(LD1) sorted marker genes 
sort(abs(lda_result$scaling[,"LD1"]))

########################################################################################
# Apply LDA on KOs:  
########################################################################################

ko_names <- c("EED KO E8.5", "RNF2 KO E8.5", "KDM6A KO E8.5", "KDM6A/B UTY TKO E8.5") 
combined_results <- list()

for(i in 1:length(ko_objs)){

  test <- get(ko_objs[i])
  test_subset <- subset(test, features = selected_genes)
  name <- ko_names[i]

  # Extract data 
  test_for_lda <- t(as.matrix(GetAssayData(test_subset, slot = "data")))
  class_labels <- test_subset@meta.data$stage
  test_lda_prediction <- predict(lda_result, test_for_lda)
  str(test_lda_prediction)

  test_lda_scales <- as.data.frame(test_lda_prediction$x)
  test_lda_scales$stage <- class_labels
  test_lda_scales$stage <- factor(test_lda_scales$stage, levels = names(named_stages), labels = named_stages)

  test_subset <- AddMetaData(test_subset, metadata = as.data.frame(test_lda_scales[1:5]))
  test_subset <- AddMetaData(test_subset, metadata = as.data.frame(test_lda_prediction$posterior))

  # Compute weighted sum of stages
  weighted_sums <- apply(test_lda_prediction$posterior, 1, function(x) sum(x * stages_numeric))
  test_subset <- AddMetaData(test_subset, metadata = as.data.frame(weighted_sums))
  test_subset$name <- name
  combined_results[[i]] <- as.data.frame(test_subset@meta.data)

}

# Combine all the data frames in the list into one
classified_kos <- do.call(rbind, combined_results)

y_axis_order <- c("RNF2 KO E8.5", "EED KO E8.5", "KDM6A KO E8.5", "KDM6A/B UTY TKO E8.5")
classified_kos$name <- factor(classified_kos$name, levels = y_axis_order)

# female and male subsets
female_ko <- classified_kos[classified_kos$sex == "female", ]
male_ko <- classified_kos[classified_kos$sex == "male", ]

# Amount of cells per lineage
table(female_ko$cluster_germlayer, female_ko$name)
table(male_ko$cluster_germlayer, male_ko$name)

lineage_colors <- c(
  "Epiblast" = "#000000",
  "Blood" = "#654024",
  "Meso" = "#945BC7",
  "XEcto" = "#EE615D",
  "XEndo" = "#EC8F00",
  "XMeso" = "#68DEFF",
  "Ecto" = "#38A983",
  "Endom" = "#F7E500",
  "PGC" = "#E22334"
)

# Boxplots for single cell stage predictions
f_plot <- ggplot(female_ko, aes(x=weighted_sums, y=cluster_germlayer, color=cluster_germlayer, fill=cluster_germlayer)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA,color = "black", show.legend= FALSE) +
  geom_point(alpha=0.8, position=position_jitter(height=0.1)) + 
  labs(x="Predicted stage", y= "", color = "Lineage", title="female") +
  facet_grid(name~., switch="y") +
  scale_color_manual(values = lineage_colors) +
  scale_fill_manual(values = lineage_colors, guide="none") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20, color = "black"),  
    axis.text.y = element_blank(),  
    strip.text = element_text(size = 20),  
    axis.title = element_text(size=20),  
    title = element_text(size=20),  
    legend.text = element_text(size=20), 
    legend.title = element_text(size=20) 
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5), reverse=TRUE)  
  )
ggsave("LDA_female_KDM_core_KO_staging.pdf", f_plot, height=14, width=10)

# female as violin plot:
f_plot <- ggplot(female_ko, aes(x=weighted_sums, y=cluster_germlayer, color=cluster_germlayer, fill=cluster_germlayer)) +
  geom_violin(alpha=0.5, color = "black", show.legend= FALSE) +
  geom_point(alpha=0.8, position=position_jitter(height=0.1)) + 
  labs(x="Predicted stage", y= "", color = "Lineage", title="female") +
  facet_grid(name~., switch="y") +
  scale_color_manual(values = lineage_colors) +
  scale_fill_manual(values = lineage_colors, guide="none") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20, color = "black"),  
    axis.text.y = element_blank(),  
    strip.text = element_text(size = 20),  
    axis.title = element_text(size=20),  
    title = element_text(size=20),  
    legend.text = element_text(size=20),  
    legend.title = element_text(size=20)  
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5), reverse=TRUE)  
  )
ggsave("LDA_female_KDM_core_KO_staging_violin.pdf", f_plot, height=14, width=10)

m_plot <- ggplot(male_ko, aes(x=weighted_sums, y=cluster_germlayer, color=cluster_germlayer, fill=cluster_germlayer)) +
  geom_boxplot(alpha=0.5, outlier.shape = NA,color = "black", show.legend= FALSE) +
  geom_point(alpha=0.8, position=position_jitter(height=0.1)) + 
  labs(x="Predicted stage", y= "", color = "Lineage", title="male") +
  facet_grid(name~., switch="y") +
  scale_color_manual(values = lineage_colors) +
  scale_fill_manual(values = lineage_colors, guide="none") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(size = 20, color = "black"), 
    axis.text.y = element_blank(),  
    strip.text = element_text(size = 20),  
    axis.title = element_text(size=20),  
    title = element_text(size=20),  
    legend.text = element_text(size=20),  
    legend.title = element_text(size=20)  
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 5), reverse=TRUE)  
  )
ggsave("LDA_male_KDM_core_staging.pdf", m_plot, height=14, width=10)

# Boxplots for embryo staging 
f_box <- ggplot(female_ko, aes(x=weighted_sums, y=name))+
  geom_boxplot(fill="darksalmon", outlier.size=0.5, alpha=0.5) + 
  labs(x="Predicted stage", y= "",  title="female")+
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(size = 14, color="black"), 
    axis.text.y = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 14), axis.title = element_text(size=16), title=element_text(size=16) )
ggsave("LDA_female_KDM_core_staging_box.pdf", f_box, width=7, height=3)

m_box <- ggplot(male_ko, aes(x=weighted_sums, y=name))+
  geom_boxplot(fill="cadetblue", outlier.size=0.5, alpha=0.5) + 
  labs(x="Predicted stage", y= "",  title="male")+
  theme_classic() +
  theme(legend.position = "right", axis.text.x = element_text(size = 14, color="black"),  
    axis.text.y = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 14), axis.title = element_text(size=16), title=element_text(size=16) )
ggsave("LDA_male_KDM_core_staging_box.pdf", m_box, width=7, height=3)

f_mean <- female_ko %>% group_by(name,embryo) %>% summarise(mean = mean(weighted_sums))
m_mean <- male_ko %>% group_by(name,embryo) %>% summarise(mean = mean(weighted_sums))

m_box <- ggplot(m_mean, aes(x=mean, y=factor(name, rev(levels(factor(name))))))+
  geom_boxplot(fill="cadetblue", outlier.size=0.5, alpha=0.5) + 
  labs(x="Predicted stage", y= "",  title="male")+
  theme_classic() +
  xlim(c(6.5,9.0)) +
  theme(legend.position = "right", axis.text.x = element_text(size = 14, color="black"), 
    axis.text.y = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 14), axis.title = element_text(size=16), title=element_text(size=16) )
ggsave("LDA_male_KDM_core_staging_box_embryos.pdf", m_box, width=7, height=3)

f_box <- ggplot(f_mean, aes(x=mean, y=factor(name, rev(levels(factor(name))))))+
  geom_boxplot(fill="darksalmon", outlier.size=0.5, alpha=0.5) + 
  labs(x="Predicted stage", y= "",  title="female")+
  theme_classic() +
  xlim(c(6.5,9.0)) +
  theme(legend.position = "right", axis.text.x = element_text(size = 14, color="black"), 
    axis.text.y = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 14), axis.title = element_text(size=16), title=element_text(size=16) )
ggsave("LDA_female_KDM_core_staging_box_embryos.pdf", f_box, width=7, height=3)

# Stacked plot arranged by pseudotime (stage prediction)
# calculate cell state poportions
classified_kos$embryo_name <- paste(classified_kos$embryo, classified_kos$name)
all_combinations <- expand.grid(embryo = classified_kos$embryo_name, cluster_names = cluster_names)
total_cells <- classified_kos %>% group_by(embryo_name) %>%  summarise(total = n(), .groups = 'drop')

proportions <- classified_kos %>%
  left_join(total_cells, by = "embryo_name") %>%
  group_by(cluster_names, embryo_name, total, sex) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)
proportions_complete <- all_combinations %>%
  left_join(proportions, by = c("cluster_names", "embryo"="embryo_name")) %>%
  replace_na(list(count = 0, proportion = 0)) 
germlayer_map <- setNames(cluster_germlayer, cluster_names)
proportions_complete <- proportions_complete %>%
  mutate(cluster_germlayer = germlayer_map[cluster_names]) 
  arrange(cluster_germlayer, cluster_names) 
proportions_complete$cluster_names <- factor(proportions_complete$cluster_names, levels = unique(proportions_complete$cluster_names))
means <- classified_kos %>% group_by(name,embryo) %>% summarise(mean = mean(weighted_sums))
means$embryo_name <- paste(means$embryo, means$name, sep=" ")
merged_data <- proportions %>%
  left_join(means %>% dplyr::select(embryo_name,mean), 
            by = c("embryo_name")) %>%
  mutate(proportion = if_else(is.na(proportion), 0, proportion))
merged_data <- merged_data %>% mutate(cluster_germlayer = germlayer_map[cluster_names]) 
head(merged_data)
sum <- merged_data %>% group_by(embryo_name) %>% summarize(wm=sum(proportion))
merged_data <- merged_data %>%
  arrange(cluster_germlayer, cluster_names) 
merged_data$cluster_names <- factor(merged_data$cluster_names, levels = unique(merged_data$cluster_names))
merged_data <- merged_data %>%
  arrange(mean)
merged_data$embryo_name <- factor(merged_data$embryo_name, levels = unique(merged_data$embryo_name))

embryo_sex_label <- merged_data %>%
  group_by(embryo_name, mean, name) %>%  
  summarize(sex_label = substr(first(sex),1,1), .groups = 'drop')  
embryo_sex_label$embryo_name <- factor(embryo_sex_label$embryo_name, levels = unique(embryo_sex_label$embryo_name))

test_stacked_plot <- ggplot(merged_data, aes(x = factor(round(mean,3)), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryonic stage prediction", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~name, scales="free",space="free_x")+
        theme(legend.position = "none", axis.text.x=element_text(size=2,angle=45, hjust=1))+
        geom_text(data=embryo_sex_label,aes(x=factor(round(mean,3)), label = sex_label, y = 101),size=5,vjust = 0.25, inherit.aes = FALSE) +
        theme(axis.text.y=element_text(size = 30, color="black"),axis.text.x=element_text(size = 20, color="black"),  # Increase y axis tick text size
        strip.text = element_text(size = 30), axis.title = element_text(size=30), legend.text = element_text(size=30), legend.title = element_text(size = 30)) 
ggsave("LDA_KDM_core_stacked_cellstates_ordered_pred_stage_mean.pdf",width=21, test_stacked_plot)

