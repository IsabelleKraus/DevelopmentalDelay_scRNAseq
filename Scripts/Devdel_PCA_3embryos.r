library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(future)
library(ggExtra)
library(cowplot)
library(caret)
library(Metrics)
library(patchwork)
library(viridis)

load("WT.Robj")

# Annotations
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

stages <- unique(WT@meta.data$stage)
named_colors <- setNames(cluster_colors, cluster_names)
# Get cell names
all_cells <- colnames(WT@assays$RNA@counts)
unique_stages <- unique(WT@meta.data$stage)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, unique_stages)
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames( stage_colors, c("6.5","7.0","7.5","8.0","8.5","9.0"))

#####################################################
# Stacked cell states (ordered by embryo number)
#####################################################
metadata <- WT@meta.data 
metadata <- metadata %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))
all_combinations <- expand.grid(embryo_stage = unique(metadata$embryo_stage), 
                                cluster_names = cluster_names)

# calculate median cell state poportions 
total_cells <- metadata %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
proportions <- metadata %>%
  left_join(total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, embryo_stage, total, sex) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

proportions_complete <- all_combinations %>%
  left_join(proportions, by = c("cluster_names", "embryo_stage")) %>%
  replace_na(list(count = 0, proportion = 0))

proportions <- proportions %>% separate(embryo_stage, into = c("embryo", "stage"), sep = "__")
proportions <- proportions %>% mutate(embryo_num = as.integer(gsub("[^0-9]", "",embryo)))

germlayer_map <- setNames(cluster_germlayer, cluster_names)

# Update the dataframe by correctly assigning germlayers and order 
proportions <- proportions %>%
  mutate(cluster_germlayer = germlayer_map[cluster_names])%>%
  mutate(e_day = named_stages[stage]) %>%
  arrange(cluster_germlayer, cluster_names)

proportions$sex_label <- substr(proportions$sex, 1,1)
proportions$cluster_names <- factor(proportions$cluster_names, levels = unique(proportions$cluster_names))

embryo_sex_label <- proportions %>%
  group_by(embryo_num, e_day) %>%  # Group by both number and stage
  summarize(sex_label = first(sex_label), .groups = 'drop')  # Get the first sex label for each group

# Stacked plot
embryos_stacked_plot <- ggplot(proportions, aes(x = factor(embryo_num), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        geom_text(data= embryo_sex_label,aes(x=factor(embryo_num), label = sex_label, y = 101),vjust = 0.25, inherit.aes = FALSE, size=4) +
        facet_grid(~e_day, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme( legend.position = "none", axis.text.x = element_text(size=11,color="black",angle=45, hjust=1),  # Increase x axis tick text size
          axis.text.y = element_text(size = 20, color="black"),  # Increase y axis tick text size
          strip.text = element_text(size = 20), axis.title = element_text(size=20))
ggsave("embryos_stacked_cellstates.pdf",width=14,height=5, embryos_stacked_plot)

#################################################################################################################################################
# Select cells of one whole embryo per time point as test and the other embryos as train for the PCA staging  
# testembryos: embryo1, embryo4, embryo8
testembryos <- c("embryo1", "embryo4", "embryo6")
test <- subset(WT@meta.data, embryo == testembryos)
train <- subset(WT@meta.data, !(embryo %in% testembryos))

####################################################################################################################
# Training of PCA 
####################################################################################################################
train$embryo_stage <- paste(train$embryo, train$stage, sep = "__")

all_combinations <- expand.grid(embryo_stage = unique(train$embryo_stage), 
                                cluster_name = cluster_names)

# median cell state poportions 
total_cells <- train %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# proportions for each combination of cluster_name, stage, and embryo
proportions <- train %>%
  left_join(total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

# To get PC1 distribution of cell state proportions in each train embryo 
proportions_complete <- all_combinations %>%
  left_join(proportions, by = c("cluster_name" = "cluster_names", "embryo_stage")) %>%
  replace_na(list(count = 0, proportion = 0))

# Calculate median proportions for every stage and cluster_name over embryo: deviding s
median_proportions <- proportions_complete %>% separate(embryo_stage, into = c("embryo", "stage"), sep = "__") %>%
  group_by(stage, cluster_name) %>%
  summarise(median_proportion = median(proportion), .groups = 'drop')

sum <- median_proportions %>%
  group_by(stage) %>%
  summarise(sum = sum(median_proportion) )

# Normalize median proportions
median_proportions <- median_proportions %>%
  left_join(sum, by = "stage") %>%
  mutate(normalized_median_proportion = median_proportion / sum) %>%
  select(-sum) # Removing the sum column after normalization

# Binary Information of cell states

train_prepared <- train %>%
  distinct(embryo_stage, cluster_names) %>%
  mutate(presence = 1)

# Merge with metadata to check presence
binary_presence <- all_combinations %>%
  left_join(train_prepared, by = c("embryo_stage", "cluster_name" = "cluster_names")) %>%
  replace_na(list(presence = 0)) %>%
  separate(embryo_stage, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo, stage, cluster_name, binary_presence = presence)

# Calculate median proportions for every stage and cluster_name over embryos
median_binary <- binary_presence %>% group_by(stage, cluster_name) %>%
  summarise(median_binary = median(binary_presence), .groups = 'drop')
 
# PCA for each stage, on each embryo (median_proportion vector and on binary info vectors)

median_proportions_wide <- median_proportions %>% select(-median_proportion) %>%
  pivot_wider(names_from = cluster_name, values_from = normalized_median_proportion, values_fill = list(normalized_median_proportion = 0.0)) 
  
# Identify columns with near zero variance
nzv_prop <- nearZeroVar(median_proportions_wide)
# no near zero variance in this case 

# Remove near zero variance columns from the dataset (if there are nzvs)
median_proportions_wide_clean <- median_proportions_wide %>% select(-nzv_prop)
median_proportions_wide_clean <- median_proportions_wide_clean %>% select(-stage)

# PCA on median proportions
set.seed(42)
pca_median_proportions <- prcomp(median_proportions_wide_clean, retx=TRUE, center = TRUE, scale. = TRUE)
saveRDS(pca_median_proportions, "pca_median_proportions.rds")

pca_prop_df <- as.data.frame(pca_median_proportions$x)
pca_prop_df$stage <- stages

binary_presence_wide <- median_binary %>%
  pivot_wider(names_from = cluster_name, values_from = median_binary, values_fill = list(median_binary = 0)) %>%
  select(-stage) 

# Remove near zero variance columns from the dataset
nzv_bin <- nearZeroVar(binary_presence_wide)
# 9 22 32 34 39
binary_presence_wide_clean <- binary_presence_wide[, -nzv_bin]

set.seed(42)
pca_binary_presence <- prcomp(binary_presence_wide_clean, retx=TRUE, center = TRUE, scale. = TRUE)
saveRDS(pca_binary_presence, "pca_binary_presence.rds")

df_binary_presence <- data.frame(Stage = stages, PC1_Binary = pca_binary_presence$x[,"PC1"], PC2_Binary = pca_binary_presence$x[,"PC2"])
df_median_proportions <- data.frame(Stage = stages, PC1_Proportion = pca_median_proportions$x[,"PC1"], PC2_Proportion = pca_median_proportions$x[,"PC2"])
df_binary_presence$e_days <- e_days
df_median_proportions$e_days <- e_days

# Combine the dataframes
combined_df <- merge(df_binary_presence, df_median_proportions, by = "Stage")

# Calculate variances 
variance_bin <- pca_binary_presence$sdev^2
variance_prop <- pca_median_proportions$sdev^2
total_variance_bin <- sum(variance_bin)
total_variance_prop <- sum(variance_prop)

# Calculate proportion (percentage) of variance explained by each component 
proportion_variance_bin <- (variance_bin / total_variance_bin)*100
proportion_variance_prop <- (variance_prop / total_variance_prop)*100

unique_stages <- unique(WT@meta.data$stage)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, unique_stages)
combined_df$e_days <- e_days

# Define a set of distinct colors for each stage
color_mapping2 <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
                          e_days)

scatter <- ggplot(combined_df, aes(x=PC1_Binary, y=PC1_Proportion, label=e_days))+
  geom_point(shape=8, aes(colour=e_days))+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  theme_minimal()+
  labs(colour="Stage")+
  scale_color_manual(values=color_mapping2, )+
  theme_classic()+
  geom_text(hjust=0.5, vjust=2)
ggsave("PC1_values_median_WT_embryos_normalized.pdf", scatter)

# Scatter-Plots for 2 principal components 
scatter_bin <- ggplot(df_binary_presence, aes(x=PC1_Binary, y=PC2_Binary, label=e_days))+
  geom_point()+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC2 binary (", round(proportion_variance_bin[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  geom_text(hjust=0.5, vjust=2)
ggsave(paste0("PC1_PC2_values_median_WT_binary.pdf"), scatter_bin)

scatter_prop <- ggplot(df_median_proportions, aes(x=PC1_Proportion, y=PC2_Proportion, label=e_days))+
  geom_point()+
  xlab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  ylab(paste0("PC2 proportions (", round(proportion_variance_prop[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  geom_text(hjust=0.5, vjust=2)
ggsave(paste0("PC1_PC2_values_normalized_median_WT_proportions.pdf"), scatter_prop)


############################################################################################################################################
# Distributions of first principal components for each trainings embryo (for different stages)

binary_presence$embryo_stage <- paste(binary_presence$embryo,binary_presence$stage, sep="__")
# Reshape and perform PCA on binary presence
binary_presence_wide2 <- binary_presence %>%
  group_by(embryo_stage) %>%
  pivot_wider(names_from = cluster_name, values_from = binary_presence, values_fill = list(binary_presence = 0))

proportions_complete <- proportions_complete %>%
  mutate(
    embryo = stringr::str_extract(embryo_stage, "embryo\\d+"), # Extract 'embryo' followed by digits
    stage = stringr::str_extract(embryo_stage, "WT_\\d+") # Extract 'WT_' followed by digits
  )

proportions_wide2 <- proportions_complete %>%
  select(-total, -count) %>%
  group_by(embryo_stage) %>%
  pivot_wider(names_from = cluster_name, values_from = proportion, values_fill = list(proportion = 0)) 

embryo_stage_prop <- proportions_wide2$embryo_stage
embryo_stage_bin <- binary_presence_wide2$embryo_stage
stages_bin <- binary_presence_wide2$stage
stages_prop <- proportions_wide2$stage

proportions_wide2_clean <- proportions_wide2 %>% select(-stage, -embryo, -embryo_stage)
binary_wide2_clean <- binary_presence_wide2 %>% select(-stage, -embryo, -embryo_stage)
binary_wide2_clean <-  binary_wide2_clean %>% ungroup() %>% select(-embryo_stage)
proportions_wide2_clean <- proportions_wide2_clean %>% ungroup() %>% select(-embryo_stage)

proportions_wide2_clean <- proportions_wide2_clean %>%
  select(all_of(names(median_proportions_wide_clean)))
binary_wide2_clean <- binary_wide2_clean %>%
  select(all_of(names(binary_presence_wide_clean)))

# PCA on all train embryos 

# prediction of PCs for embryos of train WT_90 -> loadings / rotations of WT_90 (train medians) were used 
pca_train_proportions <- predict(pca_median_proportions, newdata=proportions_wide2_clean)
pca_train_binary <-  predict(pca_binary_presence, newdata=binary_wide2_clean)  
pca_prop_df2 <- as.data.frame(pca_train_proportions)
pca_bin_df2 <- as.data.frame(pca_train_binary)

# Add embryo and stage information back to it 
pca_prop_df2$embryo_stage <- embryo_stage_prop
pca_prop_df2$stage <- stages_prop
pca_bin_df2$embryo_stage <- embryo_stage_bin

train_combined_df <- 
  merge(data.frame(embryo_stage=pca_prop_df2$embryo_stage, PC1_Proportion = pca_prop_df2[,"PC1"], stage =pca_prop_df2$stage),data.frame(embryo_stage = pca_bin_df2$embryo_stage, PC1_Binary = pca_bin_df2[,"PC1"]), by="embryo_stage")%>%
  mutate(e_day = named_stages[stage]) 

# PC1 Distribution Plots: 
# Plot for PC1_Proportion
train_combined_df$day <- substr(train_combined_df$e_day, 2,4)

pc1_prop <- ggplot(train_combined_df, aes(x = day, y = PC1_Proportion)) +
  geom_boxplot() +
  theme_classic() +
  labs( x = "Embryonic stage (E)", y = "PC1", title="Proportion") +
  theme(
    axis.text.x = element_text(size = 24, color="black", margin = margin(t = 10)), 
    plot.title = element_text(size=24, margin = margin(b = 20)),  # Increase x axis tick text size
    axis.text.y = element_text(size = 24, color="black",margin = margin(r = 10)), # Increase y axis tick text size
    strip.text = element_text(size = 24), 
    axis.title.x = element_text(size=24,margin = margin(t = 10)),
    axis.title.y = element_text(size=24,margin = margin(r = 10)),
    )
# Plot for PC1_Binary
pc1_bin <- ggplot(train_combined_df, aes(x = day, y = PC1_Binary)) +
  geom_boxplot() +
  theme_classic() + 
  labs(x = "Embryonic stage (E)", y = "PC1", title="Binary (presence or absence)") +
 theme(
    axis.text.x = element_text(size = 24, color="black", margin = margin(t = 10)), 
    plot.title = element_text(size=24, margin = margin(b = 20)),  # Increase x axis tick text size
    axis.text.y = element_text(size = 24, color="black",margin = margin(r = 10)), # Increase y axis tick text size
    strip.text = element_text(size = 24), 
    axis.title.x = element_text(size=24,margin = margin(t = 10)),
    axis.title.y = element_text(size=24,margin = margin(r = 10))
    )
ggsave(paste0("PC1_distr_prop.pdf"),width=7,height=7, pc1_prop)
ggsave(paste0("PC1_distr_bin.pdf"),width=7, pc1_bin)


#############################################################################################################################

# Plot median cell state proportions (stacked barplot) and binary information for each stage 
# Add germlayer information to the dataframe

# Create a named vector for mapping cluster names to germlayers
germlayer_map <- setNames(cluster_germlayer, cluster_names)

# Update the dataframe by correctly assigning germlayers and order 
median_binary <- median_binary %>%
  mutate(cluster_germlayer = germlayer_map[cluster_name])%>%
  mutate(e_day = named_stages[stage]) %>%
  arrange(cluster_germlayer, cluster_name)
  
median_proportions <- median_proportions %>%
  mutate(cluster_germlayer = germlayer_map[cluster_name]) %>% 
  mutate(e_day = named_stages[stage]) %>%
  arrange(cluster_germlayer, cluster_name)

median_binary$cluster_name <- factor(median_binary$cluster_name, levels = unique(median_binary$cluster_name))
median_proportions$cluster_names <- factor(median_proportions$cluster_name, levels = unique(median_proportions$cluster_name))
median_proportions$day <- substr(median_proportions$e_day, 2,4)

# stacked plots with WT as reference 
median_stacked_plot <- ggplot(median_proportions, aes(x = day, y = normalized_median_proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryonic stage (E)", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", 
          axis.text.x = element_text(size = 24, color="black"), 
          axis.text.y = element_text(size = 24, color="black"), # Increase y axis tick text size
          strip.text = element_text(size = 24), 
          axis.title.x = element_text(size=24,margin = margin(t = 10)),
          axis.title.y = element_text(size=24))
ggsave(paste0("median_stacked_cellstates_normalized.pdf"),width=7, median_stacked_plot)

median_binary$combined_cluster <- paste(median_binary$cluster_name, median_binary$cluster_germlayer, sep = "_")
median_binary$combined_cluster <- factor(median_binary$combined_cluster, levels = unique(median_binary$combined_cluster))
median_binary$day <- substr(median_binary$e_day,2,4)
median_binary$cluster_num <- sub("^([0-9]+)_.*", "\\1", median_binary$cluster_name)

# binary tile plots 
binary_plot <- ggplot(median_binary, aes(x = day, y = cluster_num, fill = factor(median_binary))) +
  geom_tile(color = "grey") + # Add tiles with a grey border
  scale_fill_manual(values = c("1" = "darkred", "0" = "navy", "0.5" = "#460040"),  labels = c("0" = "0 (absent)", "1" = "1 (present)")) + # Set fill colors
  theme_minimal() + # Use a minimal theme
  facet_grid(cluster_germlayer ~ ., scales = "free_y", space = "free_y")+ # Group by germlayer with space 
  labs(x = "Embryonic stage (E)", y = "Cell State", fill = "Median Binary") + # Add labels
  theme(axis.text = element_text(size=14, color="black"), # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(), 
        strip.text.y= element_text(size=16, angle=0, hjust=0),
          axis.title.x = element_text(size=16,margin = margin(t = 10)),
          axis.title.y = element_text(size=16),# Remove minor grid lines
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))+
    guides(color=guide_legend(title="Median binary", override.aes=list(size=3)))
ggsave(paste0("median_binary_plot.pdf"),width=8, height=14, binary_plot)

##############################################################################################################################################################
# Test (Validation)
##############################################################################################################################################################

test <- test %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))

length(unique(test$embryo_stage))
length(unique(train$embryo_stage))

# Stage proportions (Cells): 
stage_proportions_test <- test %>% group_by(stage) %>% summarise(proportion = n()/dim(test)[1])
stage_proportions_test
stage_proportions_train <- train %>% group_by(stage) %>% summarise(proportion = n()/dim(train)[1])
stage_proportions_train

test_all_combinations <- expand.grid(embryo_stage = unique(test$embryo_stage), 
                                cluster_name = cluster_names)

# Calculate median cell state poportions 

test_total_cells <- test %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
test_proportions <- test %>%
  left_join(test_total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total, sex) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

test_proportions_complete <- test_all_combinations %>%
  left_join(test_proportions, by = c("cluster_name" = "cluster_names", "embryo_stage")) %>%
  replace_na(list(count = 0, proportion = 0)) %>%
  mutate(split=embryo_stage)%>%
  separate(split, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo_stage, embryo, stage, cluster_name, proportion = proportion, sex)

test_proportions_wide <- test_proportions_complete %>%
  select(embryo_stage,stage, cluster_name, proportion) %>%  # Select only the necessary columns
  distinct() %>%  # Ensure there are no duplicate rows
  pivot_wider(names_from = cluster_name, values_from = proportion, values_fill = list(proportion = 0.0))

test_embryo_stage <- test_proportions_wide$embryo_stage

# Extract column names from median_proportions_wide
column_order <- colnames(median_proportions_wide)

# Reorder columns of test_proportions_wide to match
test_proportions_wide <- test_proportions_wide %>% 
  select(all_of(column_order)) 
# remove same cluster_name as in train median proportions
#test_proportions_wide_clean <- test_proportions_wide[, -nzv_prop] 
test_proportions_wide_clean <- test_proportions_wide %>% select(-stage)

# prediction of PCs for test dataset WT_10 -> loadings / rotations of WT_90 were used 
pred_prop <- predict(pca_median_proportions, newdata=test_proportions_wide_clean)  
pred_prop <- as.data.frame(pred_prop)
# Add embryo and stage information back to it 
pred_prop$embryo_stage <- test_embryo_stage

# Stacked Plot for each test embryo: 
test_proportions_complete <- test_proportions_complete %>%
  mutate(cluster_germlayer = germlayer_map[cluster_name]) %>% 
  mutate(e_day = named_stages[stage]) %>%
  arrange(stage, cluster_germlayer, cluster_name) %>%
  mutate(order_number = (match(embryo_stage, unique(embryo_stage))))

test_proportions_complete$cluster_name <- factor(test_proportions_complete$cluster_name, levels = unique(test_proportions_complete$cluster_name))
test_proportions_complete$embryo_num <- as.integer(gsub("[^0-9]", "",test_proportions_complete$embryo))
test_proportions$embryo_num <- as.integer(gsub("[^0-9]", "",test_proportions$embryo))

embryo_sex_label <- test_proportions %>%
  group_by(embryo_num, stage) %>%  # Group by both order_number and stage
  summarize(sex_label = substr(first(sex),1,1), .groups = 'drop')  # Get the first sex label for each group
embryo_sex_label <- embryo_sex_label %>% mutate(e_day = named_stages[stage]) 

# stacked plots with WT as reference 
test_stacked_plot <- ggplot(test_proportions_complete, aes(x = factor(embryo_num), y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        geom_text(data= embryo_sex_label,aes(x=factor(embryo_num), label = sex_label, y = 101),vjust = 0.25, size=5, inherit.aes = FALSE) +
        facet_grid(~e_day, scales="free",space="free_x")+
          theme(
            legend.position = "none",
            axis.text.x = element_text(size = 16, color="black"), 
            plot.title = element_text(size=24, margin = margin(b = 20)),  # Increase x axis tick text size
            axis.text.y = element_text(size = 24, color="black"), # Increase y axis tick text size
            strip.text = element_text(size = 24), 
            axis.title.x = element_text(size=24,margin = margin(t = 10)),
            axis.title.y = element_text(size=24)
        )
ggsave(paste0("test_stacked_cellstates.pdf"),width=14, test_stacked_plot)


# Binary Information of cell states 

test_prepared <- test %>%
  distinct(embryo_stage, cluster_names) %>%
  mutate(presence = 1)

test_binary_presence <- test_all_combinations %>%
  left_join(test_prepared, by = c("embryo_stage", "cluster_name" = "cluster_names")) %>%
  replace_na(list(presence = 0)) %>%
  mutate(split = embryo_stage) %>%
  separate(split, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo_stage, embryo, stage, cluster_name, binary_presence = presence)

# Stacked Plot for each test embryo: 
test_binary_presence <- test_binary_presence %>%
  mutate(cluster_germlayer = germlayer_map[cluster_name]) %>% 
  mutate(e_day = named_stages[stage]) %>%
  arrange(cluster_germlayer, cluster_name)

test_binary_presence$cluster_name <- factor(test_binary_presence$cluster_name, levels = unique(test_binary_presence$cluster_name))
test_binary_presence <- test_binary_presence %>%
        arrange(stage, cluster_germlayer, cluster_name) %>%
        mutate(order_number = (match(embryo_stage, unique(embryo_stage)))) %>%
        mutate(embryo_num = as.integer(sub("^embryo(\\d+).*", "\\1", embryo)))

# binary tile plots 
test_binary_plot <- ggplot(test_binary_presence, aes(x = factor(embryo_num), y = sub("^([0-9]+)_.*", "\\1", cluster_name), fill = factor(binary_presence))) +
  geom_tile(color = "grey") + # Add tiles with a grey border
  scale_fill_manual(values = c("1" = "darkred", "0" = "navy", "0.5" = "#460040"),  labels = c("0" = "0 (absent)", "1" = "1 (present)")) + # Set fill colors
  theme_minimal() + # Use a minimal theme
  facet_grid(cluster_germlayer ~  e_day , scales = "free", space = "free_y")+ # Group by germlayer with space 
  labs(x = "Embryo", y = "Cell State", fill = "") + # Add labels
  theme(axis.text = element_text(size=14, color="black"), # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size=16), 
        strip.text.y= element_text(size=16, angle=0, hjust=0),
        axis.title.x = element_text(size=16,margin = margin(t = 10)),
         axis.title.y = element_text(size=16),# Remove minor grid lines
        legend.text = element_text(size=16),
        legend.title = element_text(size=16))+
    guides(color=guide_legend(title="Median binary", override.aes=list(size=3)))
ggsave(paste0("test_binary_plot.pdf"),width=10, height=14, test_binary_plot)

test_binary_wide <- test_binary_presence %>%
  select(-c(order_number, e_day, cluster_germlayer, stage, embryo)) %>%
  pivot_wider(names_from = cluster_name, values_from = binary_presence) 
test_embryo_stage_bin <- test_binary_wide$embryo_stage
test_binary_wide <- select(test_binary_wide,-embryo_stage)

column_order <- colnames(binary_presence_wide_clean)
test_binary_wide <- select(test_binary_wide,all_of(column_order)) 

# prediction of PCs for test dataset WT_10 -> loadings / rotations of WT_90 were used 
pred_bin <- predict(pca_binary_presence, newdata=test_binary_wide)  
pred_bin <- as.data.frame(pred_bin) 
pred_bin$embryo_stage <- test_embryo_stage_bin

test_combined_df <- 
  merge(data.frame(embryo_stage=pred_prop$embryo_stage, PC1_Proportion = pred_prop[,"PC1"]),data.frame(embryo_stage = pred_bin$embryo_stage, PC1_Binary = pred_bin[,"PC1"]), by="embryo_stage") %>%
  mutate(split=embryo_stage)%>%
  separate(split, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo_stage, embryo, stage, PC1_Proportion,PC1_Binary)

test_combined_df <- test_combined_df %>%
  mutate(e_days = named_stages[stage])

# Define stages and corresponding colors
colors <- c("#A1D99B", "#7BB369", "#559E3B", "#2E8920", "#17720A", "#00441B")  # Adjust these intermediate colors as desired
# Create a named vector for color mapping
color_mapping1 <- setNames(colors, e_days)
# Define a set of distinct colors for each stage
color_mapping2 <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02"),
                          e_days)

combined_df$Group <- "train (median)"
test_combined_df$Group <- "test"
train_combined_df$Group <- "train"

# scatter plots showing PC1 values (Binary+ Proportion)
scatter1 <- ggplot(combined_df, aes(x=PC1_Binary, y=PC1_Proportion,label=e_days))+
  geom_point(aes(colour=e_days, shape=Group), size=3)+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  theme_minimal()+
  theme_classic()+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  geom_text(hjust=0.5, vjust=2)+
  labs(colour="Stage") 
ggsave(paste0("PC1_values_median_WT_embryos_a.pdf"),width=8, scatter1)

# Subset both dataframes to only include common columns
common_cols <- intersect(colnames(combined_df), colnames(test_combined_df))
combined_df_common <- combined_df[, common_cols]
test_combined_df_common <- test_combined_df[, common_cols]
# Combine the dataframes
combined_full_df <- rbind(combined_df_common, test_combined_df_common)

scatter2 <- ggplot(combined_full_df, aes(x=PC1_Binary, y=PC1_Proportion,label=e_days))+
  #geom_point(aes(shape=Group, colour=named_stage), size=1) +  # Default size for all points
  geom_point(data = subset(combined_full_df, Group == "train (median)"), aes(colour=e_days,shape=Group), size=3) +
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  theme_minimal()+
  theme_classic()+
  xlim(-6.5,6.5)+
  ylim(-8, 8)+
  #scale_size_manual(values=c("train (median)"=10, "test"=1))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  geom_text(data = subset(combined_full_df, Group == "train (median)"),
            aes(label=e_days, colour=e_days), hjust=0.5, vjust=2, check_overlap = TRUE)+
  labs(colour="Stage") 
ggsave(paste0("PC1_values_median_WT_embryos_b.pdf"),width=8, scatter2)

test_scatter <- ggplot(test_combined_df, aes(x=PC1_Binary, y=PC1_Proportion, color=e_days))+
  geom_point(aes(colour=e_days,shape=Group))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  theme_minimal()+
  theme_classic()+
  labs(colour="Stage")
ggsave(paste0("PC1_values_test_WT_embryos_cm2_b.pdf"), test_scatter+
  geom_text(aes(label=e_days),hjust=0.5, vjust=2))

together <- test_scatter+ geom_point(data=combined_df, aes(x=PC1_Binary, y=PC1_Proportion, , shape=Group)) + 
  geom_text(data=combined_df, aes(label=e_days),hjust=0.5, vjust=2)
ggsave(paste0("PC1_values_test_median_WT_embryos_a.pdf"),together)

# For marginal boxplots:
train_combined_df <- train_combined_df %>% rename("e_days"="e_day")
common_cols <- intersect(colnames(combined_full_df), colnames(train_combined_df))
combined_full_df <- rbind(combined_full_df[, common_cols] , train_combined_df[, common_cols])

blue_palette <- c(
  "E6.5" = "#ADD8E6",  # Lighter blue
  "E7.0" = "#87CEEB",
  "E7.5" = "#4682B4",
  "E8.0" = "#4169E1",
  "E8.5" = "#0000FF",
  "E9.0" = "#00008B"   # Darker blue
)
red_palette <- c(
  "E6.5" = "#FFCCCC",  # Lighter red
  "E7.0" = "#FF9999",
  "E7.5" = "#FF6666",
  "E8.0" = "#FF3333",
  "E8.5" = "#FF0000",
  "E9.0" = "#CC0000"   # Darker red
)
# Convert e_days to a factor for consistent color mapping
combined_full_df$e_days <- factor(combined_full_df$e_days)

scatter3 <- ggplot(data = subset(combined_full_df, Group == "train (median)"), 
             aes(x=PC1_Binary,y=PC1_Proportion))+
  geom_point(shape=8) +
  theme_minimal()+
  theme_classic()+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  geom_text(data = subset(combined_full_df, Group == "train (median)"),
            aes(label=e_days), hjust=0.5, vjust=2,  color="black", size=3)+
  ylim(min(combined_full_df$PC1_Proportion)-0.5,max(combined_full_df$PC1_Proportion)+0.5)+
  xlim(min(combined_full_df$PC1_Binary)-0.5,max(combined_full_df$PC1_Binary)+0.5)
ggsave(paste0("PC1_values_test_median_WT_embryos_black.pdf"),scatter3)

# assign colors and shapes
combined_full_df$group_color <- ifelse(combined_full_df$Group == "test", "#F8766D", NA)  # Assign color only for test
combined_full_df$shape <- ifelse(combined_full_df$Group == "train (median)", 8, 
                                 ifelse(combined_full_df$Group == "test", 16, 16))  # Differentiate shapes

legend_plot <- ggplot(combined_full_df, aes(x=PC1_Binary, y=PC1_Proportion)) +
  geom_point(aes(color = Group, shape = Group), size = 1) +
  scale_color_manual(name = "Group", values=c("test"="red", "train"="blue", "train (median)"="black")) +
  scale_shape_manual(name= "Group", values=c("train"=16, "test"=16, "train (median)" = 8))+
  labs()+
  theme_minimal()+
  theme_classic()
l3 <- get_legend(legend_plot)

combined_full_df$fill_color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$e_days)], blue_palette[as.character(combined_full_df$e_days)])
combined_full_df$color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$e_days)], "black")

# Create the plot
scatter3 <- ggplot(combined_full_df) +
  geom_point(aes(x=PC1_Binary, y=PC1_Proportion, color=fill_color), shape=16, size=1) + # Normal dots for train
  geom_point(data = subset(combined_full_df, Group == "train (median)"), 
             aes(x=PC1_Binary, y=PC1_Proportion),color="black", shape=8) +
  geom_text(data = subset(combined_full_df, Group == "train (median)"),
            aes(x=PC1_Binary, y=PC1_Proportion, label=e_days), hjust=0.5, vjust=2, size=3) +
  scale_color_identity() +
  theme_minimal() + 
  theme_classic() + 
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)")) +
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)")) +
  ylim(min(combined_full_df$PC1_Proportion)-0.5,max(combined_full_df$PC1_Proportion)+0.5)+
  xlim(min(combined_full_df$PC1_Binary)-0.5,max(combined_full_df$PC1_Binary)+0.5)+
  labs(shape="Group", color="Stage")

l2 <- get_legend(scatter3)
scatter3 <- scatter3 + theme(legend.position="none") 
subset <- subset(combined_full_df, Group == "train" | Group == "test")            

# Marginal boxplot of x (top panel) and y (right panel)
# For PC1_Binary
xplot <- ggboxplot(combined_full_df, x = "e_days", y = "PC1_Binary", 
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), outlier.shape = 16, alpha=0.8 ,outlier.size=1,  outlier.color="fill_color",outlier.alpha=1) +
  xlab("Stage") +
  rotate() +
  rremove("xlab") +
  ylim(min(combined_full_df$PC1_Binary)-0.5,max(combined_full_df$PC1_Binary)+0.5) +
  theme(axis.text.x = element_blank()) +
  scale_fill_identity() + # Use actual color values stored in 'fill_color'
  scale_color_identity()

# For PC1_Proportion
yplot <- ggboxplot(combined_full_df, x = "e_days", y = "PC1_Proportion",
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), alpha=0.8, outlier.shape=16, outlier.size=1, outlier_color="fill_color", outlier.alpha=1) +
  rremove("legend") +
  rremove("ylab") +
  xlab("Stage") +
  ylim(min(combined_full_df$PC1_Proportion)-0.5,max(combined_full_df$PC1_Proportion)+0.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank()) +
  scale_fill_identity() + # Use actual color values stored in 'fill_color'
  scale_color_identity()

l1 <- get_legend(xplot)
xplot <- xplot+rremove("legend")

complete <- plot_grid(xplot, NULL, scatter3, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2, 3))
# Main plot with PC train, test and train (median) scatters and marginal boxplots  
ggsave(paste0("PC1_values_scatter_boxplots_complete.pdf"),complete)
ggsave(paste0("PC1_values_scatter_boxplots_legend.pdf"),ggarrange(l1,l2,l3))

# Done: PCA on test

####################################################################################################################################################################
# Classification by euclidean distance of test PC1's to train (median) PC1's 
####################################################################################################################################################################

euclidean_distance <- function(point1, point2) {
  if (length(point1) != length(point2)) {
    stop("Points must have the same dimensions")
  }
  return(sqrt(sum((point1 - point2)^2)))
}

# test_combined_df: dataframe includes predicted PCA values, embryo and actual stage  
# x and y: PCA modes as a string (x = "PC1_Binary" or "PC1_Proportion", y = "PC1_Proportion or "PC2_Binary" or "PC2_Proportion")
# combined_df: dataframe for train (median) PCA
classification <- function(test_combined_df, x, y, combined_df) {
  
# Initialize columns in test_combined_df for each stage
for (stage in unique_stages) {
  test_combined_df[[stage]] <- NA
}
for (i in 1:nrow(test_combined_df)) {
  # Get the current point from test_combined_df
  point_test <- test_combined_df[i, c(x, y)]
  # Compute distances to each point in combined_df
  for (j in 1:nrow(combined_df)) {
    point_combined <- combined_df[j, c(x, y)]
    distance <- euclidean_distance(point_test, point_combined)
    stage_name <- combined_df$Stage[j]
    test_combined_df[i, stage_name] <- distance
  }
}
# Loop through each row of test_combined_df
for (i in 1:nrow(test_combined_df)) {
  # Extract the distances for the current row, excluding non-distance columns
  distances <- test_combined_df[i, unique_stages]
  # Find the minimum distance and the corresponding stage name
  min_distance <- min(distances, na.rm = TRUE)
  closest_stage <- names(distances)[which.min(distances)]
  test_combined_df$predicted_stage[i] <- closest_stage
}
return(test_combined_df)
}

# Perform classification on test embryos 
test_combined_df <- classification(test_combined_df, "PC1_Binary", "PC1_Proportion", combined_df)

# Classification rate for stages
correct_predictions <- sum(test_combined_df$stage == test_combined_df$predicted_stage)
total_predictions <- nrow(test_combined_df)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (PC 1 of binary pressence and proportions of cell states):", classification_rate))

# Classification Statistics (confusion matrix included)
# Ensure both predicted and actual are factors with the same levels.
all_levels <- union(levels(factor(test_combined_df$predicted_stage)), levels(factor(test_combined_df$stage)))

test_combined_df$predicted_stage <- factor(test_combined_df$predicted_stage, levels = all_levels)
test_combined_df$stage <- factor(test_combined_df$stage, levels = all_levels)

# compute the confusion matrix and other test statistics
conf_mat <- confusionMatrix(test_combined_df$predicted_stage, test_combined_df$stage)
print(conf_mat)

confMat_melt <- as.data.frame(as.table(conf_mat$table))
confMat_melt$Prediction <- factor(confMat_melt$Prediction, levels = stages)
confMat_melt$Reference <- factor(confMat_melt$Reference, levels = stages)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Prediction), y =factor(Reference), fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels=named_stages)+
  scale_y_discrete(labels=named_stages)+
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Prediction", y = "Stage", fill = "Frequency")

ggsave(paste0("PC1_values_confusion.pdf"),width=10, confMat_plot)

# For root mean square error transform discrete stage into numerical values
actual_numeric <- as.double(paste0(substr(test_combined_df$stage, 4, 4),".",substr(test_combined_df$stage, 5, 5)))
predicted_numeric <- as.double(paste0(substr(test_combined_df$predicted_stage, 4, 4),".",substr(test_combined_df$predicted_stage, 5, 5) ))
test_combined_df$actual_numeric <- actual_numeric
test_combined_df$predicted_numeric <- predicted_numeric

mae <- mae(actual_numeric, predicted_numeric)
mse <- mse(actual_numeric, predicted_numeric)

print(paste0("RMSE = ", sqrt(mse)))


###################################################################
# Test stacked barplot with information of classification and delta

test_combined_df <- test_combined_df %>%
  left_join(test_proportions_complete[,c("embryo_stage", "cluster_name", "proportion", "cluster_germlayer", "embryo_num")], by="embryo_stage")

test_combined_df$delta <- test_combined_df$predicted_numeric -  test_combined_df$actual_numeric 
test_combined_df$classification <- ifelse(test_combined_df$actual_numeric==test_combined_df$predicted_numeric,"true","false")
test_combined_df$embryo_num <- paste0(test_combined_df$embryo_num,"_", test_combined_df$e_day)
# Create a factor with levels in the order of appearance to ensure the plot respects the order
test_combined_df$Embryo_Ordered <- factor(test_combined_df$embryo_num, levels = unique(test_combined_df$embryo_num))

# stacked plots for train embryos ordered by angle in PCA bin
test_stacked_plot <- ggplot(test_combined_df, aes(x = Embryo_Ordered, y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        facet_grid(~predicted_numeric, scales="free",space="free_x")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))
ggsave("test_stacked_cellstates_prediction.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
delta_df <- test_combined_df %>% group_by(embryo_stage,embryo_num,actual_numeric, predicted_numeric, delta, classification, Embryo_Ordered) %>%  summarize()
# Color Palette
myColorPalette <- colorRampPalette(c("blue","grey", "red"))
colors <- myColorPalette(11)
# colors for all possible deltas
delta_colors <- setNames(colors, seq(-2.5,2.5,0.5))

p <- ggplot(delta_df, aes(x= Embryo_Ordered , fill = as.character(delta))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Predction-Actual", values=delta_colors) + # Color transition
  labs(x = "", y = "", fill = "Predction-Actual") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lp <- get_legend(p)
ggsave("prediction_05delta.pdf",width=21,height=2, p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())
complete <- plot_grid(test_stacked_plot,p, ncol = 1, rel_heights=c(4,1),align = "h") 
ggsave(paste0("test_stacked_05delta_complete.pdf"),width=23,complete)

# To align delta and stacked plot with actual (true) heatmap 
red_palette <- c(
  "6.5" = "#FFCCCC",  # Lighter red
  "7" = "#FF9999",
  "7.5" = "#FF6666",
  "8" = "#FF3333",
  "8.5" = "#FF0000",
  "9" = "#CC0000"   # Darker red
)

a <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual Stage", values=color_mapping) +
  labs(x = "", y = "", fill = "Actual Stage") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
la <- get_legend(a)

complete <- plot_grid(test_stacked_plot,p,a, ncol = 1, rel_heights=c(4,1,1),align = "h") 
ggsave(paste0("test_stacked_05delta_actual_complete.pdf"),width=23,complete)

b <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("false"="red","true"="olivedrab")) +
  labs(x = "", y = "", fill = "Classification") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lb <- get_legend(b)

a <- a + rremove("xlab") + rremove("legend")
p <- p + rremove("xlab") + rremove("legend")
b <- b + rremove("xlab") + rremove("legend")

complete <- b/a/p/test_stacked_plot + plot_layout(heights = c(0.2,0.2,0.2,4))
ggsave(paste0("test_stacked_05delta_actual_classification_complete.pdf"),width=23,height=9,complete)
ggsave(paste0("test_stacked_05delta_actual_classification_complete_legend.pdf"),ggarrange(lb,lp,la))

# Get embryos that are false classified  
false_embryos <- delta_df %>%
  filter(classification == "false") %>%
  select(embryo_num)

##############################################################################################################################
# Classification on PC1 and PC2 of cell state proportions (wihout binary)

# df_median_proportions: PC1 and PC2 values of median (train)
# pca_prop_df2: all PC values for all train embryo cells 
# pred_prop: all PC values for test embryo cells  
# scatter_prop: scatter plot of PC1 and PC2 median (train) proportion values

df_median_proportions$Group <- "median (train)"

df_train_pca <- data.frame(PC1_Proportion = pca_prop_df2[,"PC1"], PC2_Proportion = pca_prop_df2[,"PC2"], stage = pca_prop_df2[, "stage"], embryo_stage = pca_prop_df2[, "embryo_stage"])

df_test_pca <- data.frame(PC1_Proportion = pred_prop[, "PC1"], PC2_Proportion = pred_prop[,"PC2"],  embryo_stage = pred_prop[, "embryo_stage"])
df_test_pca <- df_test_pca %>%
  separate(embryo_stage, into = c("embryo", "stage"), sep = "__")%>%
  mutate(e_day = named_stages[stage])

df_train_pca <- df_train_pca %>% mutate(e_day = named_stages[stage])  

df_train_pca$Group <- "train"
df_test_pca$Group <- "test"
df_train_pca_prop <- df_train_pca

common_cols <- intersect(colnames(df_train_pca), colnames(df_test_pca))
# Combine the dataframes
df_combi_pca <- rbind( df_train_pca[, common_cols], df_test_pca[, common_cols])

df_medians <- df_median_proportions 
df_medians <- df_medians %>% rename("e_day"="e_days")

common_cols <- intersect(colnames(df_medians), colnames(df_combi_pca))
# Combine the dataframes
combined_full_df <- rbind( df_medians[, common_cols], df_combi_pca[, common_cols])

scatter_prop <- scatter_prop +  coord_cartesian(ylim = c(min(combined_full_df$PC2_Proportion)-0.5, max(combined_full_df$PC2_Proportion)+0.5), xlim = c(min(combined_full_df$PC1_Proportion)-0.5, max(combined_full_df$PC1_Proportion)+0.5)) 

test_scatter <- ggplot(df_test_pca, aes(x=PC1_Proportion, y=PC2_Proportion, color=e_day))+
  geom_point(aes(colour=e_day,shape=Group))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  xlab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  ylab(paste0("PC2 proportions (", round(proportion_variance_prop[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  labs(colour="Stage")
ggsave(paste0("PC1_2_values_test_WT_embryos.pdf"), test_scatter)

colnames(df_median_proportions)[colnames(df_median_proportions) == "e_days"] <- "e_day"
together <- test_scatter + 
  geom_point(data=df_median_proportions, aes(x=PC1_Proportion, y=PC2_Proportion, shape=Group)) + 
  geom_text(data=df_median_proportions, aes(label=e_day), hjust=0.5, vjust=2)
ggsave(paste0("PC1_2_test_median_WT_embryos_l.pdf"),together)

# assign colors
red_palette <- c(
  "6.5" = "#FFCCCC",  # Lighter red
  "7.0" = "#FF9999",
  "7.5" = "#FF6666",
  "8.0" = "#FF4C4C",  
  "8.5" = "#FF1919",  
  "9.0" = "#CC0000"   # Darker red
)
blue_palette <- c(
  "6.5" = "#B0E0E6",  
  "7.0" = "#87CEEB",
  "7.5" = "#5F9EA0",  
  "8.0" = "#4682B4",
  "8.5" = "#1E90FF",  
  "9.0" = "#00008B"   # Darker blue
)

combined_full_df$day <- substr(combined_full_df$e_day, 2,4)
combined_full_df$fill_color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$day)], blue_palette[as.character(combined_full_df$day)])
combined_full_df$color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$day)], "black")

# Create the plot
scatter3 <- ggplot(combined_full_df) +
  geom_point(aes(x=PC1_Proportion, y=PC2_Proportion, color=fill_color),alpha=0.8, shape=16, size=3) + # Normal dots for train
  geom_point(data = subset(combined_full_df, Group == "median (train)"), 
             aes(x=PC1_Proportion, y=PC2_Proportion),color="black", shape=8, size=3) +
  geom_text(data = subset(combined_full_df, Group == "median (train)"),
            aes(x=PC1_Proportion, y=PC2_Proportion, label=day), hjust=0.5, vjust=2, size=5) +
  scale_color_identity() +
  theme_minimal() + 
  theme_classic() + 
  xlab(paste0("PC1 (", round(proportion_variance_prop[1]), "%)")) +
  ylab(paste0("PC2 (", round(proportion_variance_prop[2]), "%)")) +
  ylim(c(min(combined_full_df$PC2_Proportion)-0.5, max(combined_full_df$PC2_Proportion)+0.5)) +
  xlim(c(min(combined_full_df$PC1_Proportion)-0.5, max(combined_full_df$PC1_Proportion)+0.5))+
  labs(shape="Group", color="Embryonic stage (E)")+
  theme(axis.text = element_text(size = 14, color="black", margin = margin(t = 10)),
     axis.title = element_text(size=14,margin = margin(t = 10)), 
     legend.position="right",
     legend.title = element_text(size=14),
     legend.text = element_text(size=14)) 

scatter3 <- scatter3 + theme(legend.position="none") 

subset <- subset(combined_full_df, Group == "train" | Group == "test")            
# Marginal boxplot of x (top panel) and y (right panel)
#xplot <- ggMarginal(scatter3,data = subset, x = "e_days", y = "PC1_Proportion",type="boxplot", color=Group, fill=Group, margins="x")

# Marginal boxplot of x (top panel) and y (right panel)
# For PC1_Proportion
xplot <- ggboxplot(combined_full_df, x = "day", y = "PC1_Proportion", 
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), outlier.shape = 16, alpha=0.8 ,outlier.size=1,  outlier.color="fill_color",outlier.alpha=1) +
  xlab("Embryonic stage (E)") +
  rotate() +
  rremove("xlab") +
  ylim(c(min(combined_full_df$PC1_Proportion)-0.5, max(combined_full_df$PC1_Proportion)+0.5))+
  theme(
  axis.text.x = element_blank(), 
  axis.text.y = element_text(size = 14, color="black"), # Increase y axis tick text size
  axis.title.y = element_text(size=14,margin = margin(r = 10)))+
  scale_fill_identity() + # Use actual color values stored in 'fill_colo
  scale_fill_identity() + # Use actual color values stored in 'fill_color'
  scale_color_identity()

# For PC2_Proportion
yplot <- ggboxplot(combined_full_df, x = "day", y = "PC2_Proportion",
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), alpha=0.8, outlier.shape=16, outlier.size=1, outlier_color="fill_color", outlier.alpha=1) +
  rremove("legend") +
  rremove("ylab") +
  xlab("Embryonic stage (E)") +
  ylim(c(min(combined_full_df$PC2_Proportion)-0.5, max(combined_full_df$PC2_Proportion)+0.5)) +
  theme(
  axis.text.x = element_text(angle=45, hjust=1,size = 14, color="black"), 
  axis.text.y = element_blank(), # Increase y axis tick text size
  axis.title.x = element_text(size=14,margin = margin(t = 10)))+
  scale_fill_identity() + # Use actual color values stored in 'fill_color'
  scale_color_identity()

# Adjusting specific margins for each plot
xplot <- xplot + theme(plot.margin = unit(c(2, 2, 0, 2), "pt"))  # Decrease bottom margin
yplot <- yplot + theme(plot.margin = unit(c(2, 0, 2, 2), "pt"))  # Decrease left margin
scatter3 <- scatter3 + theme(plot.margin = unit(c(0, 2, 2, 0), "pt"))  # Decrease top and right margins

complete <- plot_grid(xplot, NULL, scatter3, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2,3))
ggsave(paste0("PC1_PC2_proportion_values_scatter_boxplots_complete.pdf"),complete)

# Plot Legend:
# Create data frames for each palette
red_data <- data.frame(
  Stage = factor(names(red_palette), levels = names(red_palette)),
  Color = red_palette
)
blue_data <- data.frame(
  Stage = factor(names(blue_palette), levels = names(blue_palette)),
  Color = blue_palette
)
median_data <- data.frame(c=factor("median (train)"), d=1)

# Create dummy plots for each palette
red_plot <- ggplot(red_data, aes(x=Stage, y=1, color=Stage)) +
  geom_point(size=5) +
  scale_color_manual(values = red_palette, name = "test") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
blue_plot <- ggplot(blue_data, aes(x=Stage, y=1, color=Stage)) +
  geom_point(size=5) +
  scale_color_manual(values = blue_palette, name = "train") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
median_plot <- ggplot(median_data,aes( x=c, y=d, shape=c, color=c)) +
  geom_point(size=5) +
  theme_void() +
  scale_shape_manual(values = c("median (train)" = 8)) +
  scale_color_manual(values = c("median (train)" = "black")) +
  guides(shape = guide_legend(title = ""), color = guide_legend(title = "")) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
# Extract the legends
red_legend <- get_legend(red_plot)
blue_legend <- get_legend(blue_plot)
median_legend <- get_legend(median_plot)
ggsave(paste0("PC1_PC2_scatter_boxplots_complete_test_legend.pdf"),red_legend)
ggsave(paste0("PC1_PC2_scatter_boxplots_complete_train_legend.pdf"),blue_legend)
ggsave(paste0("PC1_PC2_scatter_boxplots_complete_median_legend.pdf"),median_legend)

############################################################################################################################################
# Using angles in PCA prop to order embryos 

# use E6.5 train median PC1 and PC2 vector as a reference for angle calculation 
ref_vec <- c(df_median_proportions[1,"PC2_Proportion"],df_median_proportions[1,"PC1_Proportion"])
df_train_pca$Angle = (atan2(df_train_pca$PC2_Proportion, df_train_pca$PC1_Proportion) * (180 / pi) + 360) %% 360
# Calculate the angle of the reference vector and normalize angles to 0-360 range
ref_angle <- (atan2(ref_vec[1], ref_vec[2]) * (180 / pi) + 360) %% 360
# Calculate relative angles
df_train_pca$Angle <- (df_train_pca$Angle - ref_angle + 360) %% 360

angle_df <- merge(proportions, df_train_pca, by="embryo_stage") %>%
  mutate(cluster_germlayer = germlayer_map[cluster_names])%>%
  arrange(cluster_germlayer, cluster_names)

angle_df$cluster_names <- factor(angle_df$cluster_names, levels = unique(angle_df$cluster_names))
angle_df <- angle_df[order(angle_df$Angle),]
angle_df$embryo_num <- paste0(as.character(parse_number(angle_df$embryo)),"_", angle_df$e_day)
angle_df$Embryo_Ordered <- factor(angle_df$embryo_num, levels = unique(angle_df$embryo_num))

# stacked plots for train embryos ordered by angle in PCA bin
angle_stacked_plot <- ggplot(angle_df, aes(x = Embryo_Ordered, y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryos ordered by angle in PCA", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete(limits = levels(angle_df$Embryo_Ordered))+
        theme( legend.position = "none", axis.text.x = element_text(size=11,color="black",angle=45, hjust=1),  
          axis.text.y = element_text(size = 20, color="black"),  
          strip.text = element_text(size = 20), 
          axis.title.x = element_text(size=20, margin =margin(t=10)),
          axis.title.y = element_text(size=20))
ggsave(paste0("angle_order_propPCA_stacked_cellstates_t.pdf"), width=14, angle_stacked_plot)

#################################################################################################################################################

# classification by euclidean distance of PC1 & PC1 prop to median (train) PC1 & PC1  
df_test_pca <- classification(df_test_pca, "PC1_Proportion", "PC2_Proportion", df_median_proportions)
df_test_pca$predicted_stage <- named_stages[df_test_pca$predicted_stage]

# Classification rate for stages
correct_predictions <- sum(df_test_pca$e_day == df_test_pca$predicted_stage)
total_predictions <- nrow(df_test_pca)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (PC 1 & 2 of cell states proportions):", classification_rate))

# Classification Statistics (confusion matrix included)
all_levels <- union(levels(factor(df_test_pca$predicted_stage)), levels(factor(df_test_pca$e_day)))

df_test_pca$predicted_stage <- factor(df_test_pca$predicted_stage, levels = all_levels)
df_test_pca$e_day <- factor(df_test_pca$e_day, levels = all_levels)
# compute the confusion matrix.
conf_mat <- confusionMatrix(df_test_pca$predicted_stage, df_test_pca$e_day)
print(conf_mat)

# Confusion Matrix plot:
confMat_melt <- as.data.frame(as.table(conf_mat$table))
confMat_melt$Prediction <- factor(confMat_melt$Prediction, levels = e_days)
confMat_melt$Reference <- factor(confMat_melt$Reference, levels = e_days)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]
confMat_melt$Prediction <- substr(confMat_melt$Prediction,2,4)
confMat_melt$Reference <- substr(confMat_melt$Reference, 2,4)

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Reference), y =factor(Prediction, levels = rev(levels(factor(Prediction)))), fill = Freq)) +
  geom_tile(color="grey") +
  scale_fill_gradient(low = "white", high = "navy") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Embryonic stage (E)", y = "Prediction (E)", fill = "Frequency")+
  theme(axis.text=element_text(size = 24, color="black"),  # Increase y axis tick text size
        strip.text = element_text(size = 24), axis.title.x = element_text(size=24,margin = margin(t = 10)),
        axis.title.y = element_text(size=24,margin = margin(r = 10)),
        legend.text = element_text(size=24), legend.title = element_text(size = 24, angle=90)) 

ggsave(paste0("PC1_PC2_proportion_confusion.pdf"),width=8, confMat_plot)

actual_numeric <- as.double(substr(df_test_pca$e_day, 2, 4))
predicted_numeric <- as.double(substr(df_test_pca$predicted_stage, 2, 4))
df_test_pca$actual_numeric <- actual_numeric
df_test_pca$predicted_numeric <- predicted_numeric

mse <- mse(actual_numeric, predicted_numeric)
print(paste0("RMSE = ", sqrt(mse)))

#################################################################
# Test stacked barplot with information of classification and delta

df_test_pca$embryo_stage <- paste0(df_test_pca$embryo, "__", df_test_pca$stage)
df_test_pca <- df_test_pca %>%
  left_join(test_proportions_complete[,c("embryo_stage", "cluster_name", "proportion", "cluster_germlayer", "embryo_num")], by="embryo_stage")
df_test_pca$delta <- df_test_pca$predicted_numeric -  df_test_pca$actual_numeric 
df_test_pca$classification <- ifelse(df_test_pca$actual_numeric==df_test_pca$predicted_numeric,"true","false")
df_test_pca$embryo_num <- paste0(df_test_pca$embryo_num,"_", substr(df_test_pca$e_day,2,4))
df_test_pca$Embryo_Ordered <- factor(df_test_pca$embryo_num, levels = unique(df_test_pca$embryo_num))

# stacked plots for train embryos ordered by angle in PCA bin
test_stacked_plot <- ggplot(df_test_pca, aes(x = Embryo_Ordered, y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        facet_grid(~predicted_numeric, scales="free",space="free_x")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))+
        theme(axis.text=element_text(size = 30, color="black"),  
        strip.text = element_text(size = 30), axis.title = element_text(size=30), legend.text = element_text(size=30), legend.title = element_text(size = 30)) 
ggsave("test_stacked_cellstates_prediction_PC12prop_t.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
delta_df <- df_test_pca %>% group_by(embryo_stage,embryo_num,actual_numeric, predicted_numeric, delta, classification, Embryo_Ordered) %>%  summarize()

# Color Palette
myColorPalette <- colorRampPalette(c("blue","grey", "red"))
# Generate a vector of colors from the palette
colors <- myColorPalette(11)
# colors for all possible deltas
delta_colors <- setNames(colors, seq(-2.5,2.5,0.5))
delta_df$predicted_numeric <- format(delta_df$predicted_numeric, nsmall = 1)
delta_df$actual_numeric <- format(delta_df$actual_numeric, nsmall = 1)

p <- ggplot(delta_df, aes(x= Embryo_Ordered , fill = as.character(delta))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Prediction error\n (Predction - Actual stage)", values=delta_colors) + # Color transition
  labs(x = "", y = "", fill = "Prediction error\n (Predction - Actual stage)") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), legend.text = element_text(size=24), legend.title = element_text(size = 24),
      axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(color = guide_legend(title = "Prediction error\n (Predction - Actual stage)",override.aes = list(size = 3))) 
lp <- get_legend(p)
ggsave("prediction_05delta_PC12prop.pdf",width=21,height=2, p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())
test_stacked_plot <- test_stacked_plot + labs(x="Embryo",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.title.x=element_text(size=30), strip.text.x = element_blank())

complete <- plot_grid(test_stacked_plot,p, ncol = 1, rel_heights=c(4,1),align = "h") 
ggsave(paste0("test_stacked_05delta_complete_PC12prop.pdf"),width=23,complete)

# To align delta and stacked plot with actual (true) heatmap 

a <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual stage (E)", values=color_mapping) +
  labs(x = "", y = "", fill = "Actual stage (E)") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
    theme(axis.text.y=element_blank(), legend.title = element_text(size = 24),  legend.text = element_text(size=24),
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm"))+
  guides(color = guide_legend(title = "Actual stage (E))",override.aes = list(size = 3))) 
la <- get_legend(a)

complete <- plot_grid(test_stacked_plot,p,a, ncol = 1, rel_heights=c(4,1,1),align = "h") 
ggsave(paste0("test_stacked_05delta_actual_complete_PC12prop.pdf"),width=23,complete)

b <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("false"="red","true"="olivedrab")) +
  labs(x = "", y = "", fill = "Classification") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  guides(color = guide_legend(title = "Classification",override.aes = list(size = 5))) +
  theme(axis.text.y=element_blank(),  legend.text = element_text(size=24),legend.title = element_text(size = 24),  strip.text = element_text(size = 30),
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lb <- get_legend(b)

a <- a + rremove("xlab") + rremove("legend")
p <- p + rremove("xlab") + rremove("legend")
b <- b + rremove("xlab") + rremove("legend")

complete <- b/a/p/test_stacked_plot + plot_layout(heights = c(0.2,0.2,0.2,4))
ggsave(paste0("test_stacked_05delta_actual_classification_complete_PC12prop.pdf"),width=23,height=9,complete)
ggsave(paste0("test_stacked_05delta_actual_classification_complete_PC12pro_legend_t.pdf"),ggarrange(lb,lp,la))

# Get embryos that are false classified  
false_embryos <- delta_df %>%
  filter(classification == "false") %>%
  select(embryo_num)
#################################################################################################################################################
# Classification using PC1 and 2 of Binary Presence 

# pca_binary_presence
df_binary_presence$Group <- "median (train)"

df_train_pca <- data.frame(PC1_Binary = pca_bin_df2[,"PC1"], PC2_Binary = pca_bin_df2[,"PC2"],  embryo_stage = pca_bin_df2[, "embryo_stage"])
df_test_pca <- data.frame(PC1_Binary = pred_bin[, "PC1"], PC2_Binary = pred_bin[,"PC2"],  embryo_stage = pred_bin[, "embryo_stage"])

df_test_pca <- df_test_pca %>%
  separate(embryo_stage, into = c("embryo", "stage"), sep = "__")%>%
  mutate(e_day = named_stages[stage])
df_train_pca <- df_train_pca %>%
  separate(embryo_stage, into = c("embryo", "stage"), sep = "__")%>%
  mutate(e_day = named_stages[stage])

df_train_pca$Group <- "train"
df_test_pca$Group <- "test"

common_cols <- intersect(colnames(df_train_pca), colnames(df_test_pca))
df_combi_pca <- rbind( df_train_pca[, common_cols], df_test_pca[, common_cols])
scatter_bin <- scatter_bin +  coord_cartesian(ylim = c(-5.5, 5.5), xlim = c(-7, 7)) 

# Marginal boxplot of x (top panel) and y (right panel)
xplot <- ggboxplot(df_combi_pca, x = "e_day", y = "PC1_Binary", 
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  xlab("Stage")+
  rotate()+
  rremove("xlab")+
  ylim(-7,7)

yplot <- ggboxplot(df_combi_pca, x = "e_day", y = "PC2_Binary",
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  rremove("legend")+
  rremove("ylab")+
  xlab("Stage")+
  ylim(-5.5,5.5)

l1 <- get_legend(xplot)
xplot <- xplot+rremove("legend")

complete <- plot_grid(xplot, NULL, scatter_bin, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2, 3))

test_scatter <- ggplot(df_test_pca, aes(x=PC1_Binary, y=PC2_Binary, color=e_day))+
  geom_point(aes(colour=e_day,shape=Group))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC2 binary (", round(proportion_variance_bin[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  labs(colour="Stage")

colnames(df_binary_presence)[colnames(df_binary_presence) == "e_days"] <- "e_day"

together <- test_scatter + 
  geom_point(data=df_binary_presence, aes(x=PC1_Binary, y=PC2_Binary, shape=Group)) + 
  geom_text(data=df_binary_presence, aes(label=e_day), hjust=0.5, vjust=2)

df_medians <- df_binary_presence 

common_cols <- intersect(colnames(df_medians), colnames(df_combi_pca))
combined_full_df <- rbind( df_medians[, common_cols], df_combi_pca[, common_cols])
combined_full_df$day <- substr(combined_full_df$e_day,2,4)
combined_full_df$fill_color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$day)], blue_palette[as.character(combined_full_df$day)])
combined_full_df$color <- ifelse(combined_full_df$Group == "test", red_palette[as.character(combined_full_df$day)], "black")

scatter3 <- ggplot(combined_full_df) +
  geom_point(aes(x=PC1_Binary, y=PC2_Binary, color=fill_color),alpha=0.8, shape=16, size=3) + # Normal dots for train
  geom_point(data = subset(combined_full_df, Group == "median (train)"), 
             aes(x=PC1_Binary, y=PC2_Binary),color="black", shape=8, size=3) +
  geom_text(data = subset(combined_full_df, Group == "median (train)"),
            aes(x=PC1_Binary, y=PC2_Binary, label=day, vjust=ifelse(day == "9.0", -0.5, 2)),
          hjust=0.5, size=5) +
  scale_color_identity() +
  theme_minimal() + 
  theme_classic() + 
  xlab(paste0("PC1 (", round(proportion_variance_bin[1]), "%)")) +
  ylab(paste0("PC2 (", round(proportion_variance_bin[2]), "%)")) +
  ylim(c(min(combined_full_df$PC2_Binary)-0.5, max(combined_full_df$PC2_Binary)+0.5)) +
  xlim(c(min(combined_full_df$PC1_Binary)-0.5, max(combined_full_df$PC1_Binary)+0.5))+
  labs(shape="Group", color="Embryonic stage (E)")+
  theme(axis.text = element_text(size = 14, color="black", margin = margin(t = 10)),
     axis.title = element_text(size=14,margin = margin(t = 10)), 
     legend.title = element_text(size=14),
     legend.text = element_text(size=14)) 

l2 <- get_legend(scatter3)
scatter3 <- scatter3 + theme(legend.position="none") 
subset <- subset(combined_full_df, Group == "train" | Group == "test")            

# Marginal boxplot of x (top panel) and y (right panel)
# For PC1_Binary
xplot <- ggboxplot(combined_full_df, x = "day", y = "PC1_Binary", 
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), outlier.shape = 16, alpha=0.8 ,outlier.size=1,  outlier.color="fill_color",outlier.alpha=1) +
  xlab("Embryonic stage (E)") +
  rotate() +
  rremove("xlab") +
  ylim(min(combined_full_df$PC1_Binary)-0.5, max(combined_full_df$PC1_Binary)+0.5) +
  theme(
  axis.text.x = element_blank(), 
  axis.text.y = element_text(size = 14, color="black"), # Increase y axis tick text size
  axis.title.y = element_text(size=14,margin = margin(r = 10)))+
  scale_fill_identity() + 
  scale_fill_identity() + 
  scale_color_identity()

# For PC2_Binary
yplot <- ggboxplot(combined_full_df, x = "day", y = "PC2_Binary",
                   color = "fill_color", fill = "fill_color", 
                   ggtheme = theme_classic(), alpha=0.8, outlier.shape=16, outlier.size=1, outlier_color="fill_color", outlier.alpha=1) +
  rremove("legend") +
  rremove("ylab") +
  xlab("Embryonic stage (E)") +
  ylim(min(combined_full_df$PC2_Binary)-0.5, max(combined_full_df$PC2_Binary)+0.5) +
  theme(
  axis.text.x = element_text(angle=45, hjust=1,size = 14, color="black"), 
  axis.text.y = element_blank(), 
  axis.title.x = element_text(size=14,margin = margin(t = 10)))+
  scale_fill_identity() + 
  scale_color_identity()

l1 <- get_legend(xplot)
xplot <- xplot+rremove("legend")

# Adjusting specific margins for each plot
xplot <- xplot + theme(plot.margin = unit(c(2, 2, 0, 2), "pt"))  # Decrease bottom margin
yplot <- yplot + theme(plot.margin = unit(c(2, 0, 2, 2), "pt"))  # Decrease left margin
scatter3 <- scatter3 + theme(plot.margin = unit(c(0, 2, 2, 0), "pt"))  # Decrease top and right margins

complete <- plot_grid(xplot, NULL, scatter3, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2, 3))

ggsave(paste0("PC1_PC2_Binary_values_scatter_boxplots_complete.pdf"),complete)


# classification by euclidean distance of PC1 & PC2 to median (train) PC1 & PC2  
df_test_pca <- classification(df_test_pca, "PC1_Binary", "PC2_Binary", df_binary_presence)
df_test_pca$predicted_stage <- named_stages[df_test_pca$predicted_stage ]

# Classification rate for stages
correct_predictions <- sum(df_test_pca$e_day == df_test_pca$predicted_stage)
total_predictions <- nrow(df_test_pca)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (PC 1 & 2 of cell states binary):", classification_rate))

# Classification Statistics (confusion matrix included)
all_levels <- union(levels(factor(df_test_pca$predicted_stage)), levels(factor(df_test_pca$e_day)))
df_test_pca$predicted_stage <- factor(df_test_pca$predicted_stage, levels = all_levels)
df_test_pca$e_day <- factor(df_test_pca$e_day, levels = all_levels)
conf_mat <- confusionMatrix(df_test_pca$predicted_stage, df_test_pca$e_day)
print(conf_mat)

# Confusion Matrix plot:
confMat_melt <- as.data.frame(as.table(conf_mat$table))
confMat_melt$Prediction <- factor(confMat_melt$Prediction, levels = e_days)
confMat_melt$Reference <- factor(confMat_melt$Reference, levels = e_days)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]
confMat_melt$Prediction <- substr(confMat_melt$Prediction, 2, 4)
confMat_melt$Reference <- substr(confMat_melt$Reference, 2, 4)

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Reference), y =factor(Prediction, levels = rev(levels(factor(Prediction)))), fill = Freq)) +
  geom_tile(color="grey") +
  scale_fill_gradient(low = "white", high = "navy") +
  theme_minimal() +
  labs(x = "Embryonic stage (E)", y = "Prediction (E)", fill = "Frequency")+
  theme(axis.text=element_text(size = 24, color="black"), 
        strip.text = element_text(size = 24), axis.title.x = element_text(size=24,margin = margin(t = 10)),
        axis.title.y = element_text(size=24,margin = margin(r = 10)),
        legend.text = element_text(size=24), legend.title = element_text(size = 24, angle=90)) 
ggsave(paste0("PC1_PC2_binary_confusion.pdf"),width=9, confMat_plot)

actual_numeric <- as.double(substr(df_test_pca$e_day, 2, 4))
predicted_numeric <- as.double(substr(df_test_pca$predicted_stage, 2, 4))
df_test_pca$actual_numeric <- actual_numeric
df_test_pca$predicted_numeric <- predicted_numeric

mae <- mae(actual_numeric, predicted_numeric)
mse <- mse(actual_numeric, predicted_numeric)
print(paste0("RMSE = ", sqrt(mse)))

###############################################################################################################################################################
# Pseudotime:
# Order Embryos based on Angle of PC vector to origin  
###############################################################################################################################################################
# Using angles in PCA bin to order embryos 

# use E6.5 train median PC1 and PC2 vector as a reference for angle calculation 
ref_vec <- c(df_binary_presence[1,"PC2_Binary"],df_binary_presence[1,"PC1_Binary"])

df_train_pca$Angle = (atan2(df_train_pca$PC2_Binary, df_train_pca$PC1_Binary) * (180 / pi) + 360) %% 360

# Calculate the angle of the reference vector and normalize angles to 0-360 range
ref_angle <- (atan2(ref_vec[1], ref_vec[2]) * (180 / pi) + 360) %% 360

# Calculate relative angles
df_train_pca$Angle <- (df_train_pca$Angle - ref_angle + 360) %% 360
df_train_pca$embryo_stage <- paste0(df_train_pca$embryo, "__", df_train_pca$stage)

angle_df <- merge(proportions, df_train_pca, by="embryo_stage")
angle_df <- angle_df %>%
  mutate(cluster_germlayer = germlayer_map[cluster_names])%>%
  arrange(cluster_germlayer, cluster_names)

angle_df$cluster_names <- factor(angle_df$cluster_names, levels = unique(angle_df$cluster_names))
angle_df <- angle_df[order(angle_df$Angle),]
angle_df$embryo_num <- paste0(as.character(parse_number(angle_df$embryo.x)),"_", substr(angle_df$e_day,2,4))
angle_df$Embryo_Ordered <- factor(angle_df$embryo_num, levels = unique(angle_df$embryo_num))

# stacked plots for train embryos ordered by angle in PCA bin
angle_stacked_plot <- ggplot(angle_df, aes(x = Embryo_Ordered, y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo (ordered by PCA)", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete(limits = levels(angle_df$Embryo_Ordered))+
        #geom_text(data= embryo_sex_label,aes(x=factor(embryo_num), label = sex_label, y = 101),vjust = 0.25, inherit.aes = FALSE) +
        #facet_grid(~e_day, scales="free",space="free_x")+
        #ggtitle(paste0("KO ", ko_name))+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))+
        theme(
          axis.text.x = element_text(size = 10, color="black", margin = margin(t = 10)), 
          plot.title = element_text(size=20, margin = margin(b = 20)),  # Increase x axis tick text size
          axis.text.y = element_text(size = 20, color="black",margin = margin(r = 10)), # Increase y axis tick text size
          strip.text = element_text(size = 20), 
          axis.title.x = element_text(size=20,margin = margin(t = 10)),
          axis.title.y = element_text(size=20,margin = margin(r = 10)),
          )
ggsave(paste0("angle_order_binaryPCA_stacked_cellstates_t1.pdf"), height=6,width=14, angle_stacked_plot)


###############################################################################################################################################################

# Test stacked barplot with information of classification and delta
###############################################################################################################################################################

df_test_pca$embryo_stage <- paste0(df_test_pca$embryo, "__", df_test_pca$stage)
df_test_pca <- df_test_pca %>%
  left_join(test_proportions_complete[,c("embryo_stage", "cluster_name", "proportion", "cluster_germlayer", "embryo_num")], by="embryo_stage")
df_test_pca$delta <- df_test_pca$predicted_numeric -  df_test_pca$actual_numeric 
df_test_pca$classification <- ifelse(df_test_pca$actual_numeric==df_test_pca$predicted_numeric,"true","false")
df_test_pca$embryo_num <- paste0(df_test_pca$embryo_num,"_", df_test_pca$e_day)
df_test_pca$Embryo_Ordered <- factor(df_test_pca$embryo_num, levels = unique(df_test_pca$embryo_num))

# stacked plots for train embryos ordered by angle in PCA bin
test_stacked_plot <- ggplot(df_test_pca, aes(x = Embryo_Ordered, y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        facet_grid(~predicted_numeric, scales="free",space="free_x")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))
ggsave("test_stacked_cellstates_prediction_PC12bin.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
delta_df <- df_test_pca %>% group_by(embryo_stage,embryo_num,actual_numeric, predicted_numeric, delta, classification, Embryo_Ordered) %>%  summarize()
# Color Palette
myColorPalette <- colorRampPalette(c("blue","grey", "red"))
colors <- myColorPalette(11)
# colors for all possible deltas
delta_colors <- setNames(colors, seq(-2.5,2.5,0.5))

p <- ggplot(delta_df, aes(x= Embryo_Ordered , fill = as.character(delta))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Predction-Actual", values=delta_colors) + # Color transition
  labs(x = "", y = "", fill = "Predction-Actual") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lp <- get_legend(p)
ggsave("prediction_05delta_PC12bin.pdf",width=21,height=2, p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())

complete <- plot_grid(test_stacked_plot,p, ncol = 1, rel_heights=c(4,1),align = "h") 
ggsave(paste0("test_stacked_05delta_complete_PC12bin.pdf"),width=23,complete)

# To align delta and stacked plot with actual (true) heatmap 

red_palette <- c(
  "6.5" = "#FFCCCC",  # Lighter red
  "7" = "#FF9999",
  "7.5" = "#FF6666",
  "8" = "#FF3333",
  "8.5" = "#FF0000",
  "9" = "#CC0000"   # Darker red
)
a <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual Stage", values=red_palette) +
  labs(x = "", y = "", fill = "Actual Stage") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
la <- get_legend(a)

complete <- plot_grid(test_stacked_plot,p,a, ncol = 1, rel_heights=c(4,1,1),align = "h") 
ggsave(paste0("test_stacked_05delta_actual_complete_PC12bin.pdf"),width=23,complete)

b <- ggplot(delta_df, aes(x = Embryo_Ordered, fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("false"="red","true"="olivedrab")) +
  labs(x = "", y = "", fill = "Classification") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lb <- get_legend(b)

a <- a + rremove("xlab") + rremove("legend")
p <- p + rremove("xlab") + rremove("legend")
b <- b + rremove("xlab") + rremove("legend")

complete <- b/a/p/test_stacked_plot + plot_layout(heights = c(0.2,0.2,0.2,4))
ggsave(paste0("test_stacked_05delta_actual_classification_complete_PC12bin.pdf"),width=23,height=9,complete)
ggsave(paste0("test_stacked_05delta_actual_classification_complete_PC12bin_legend.pdf"),ggarrange(lb,lp,la))

# Get embryos that are false classified  
false_embryos <- delta_df %>%
  filter(classification == "false") %>%
  select(embryo_num)
