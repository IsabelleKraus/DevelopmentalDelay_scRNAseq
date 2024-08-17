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
library(ggfortify)


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

# Determine the size of subsets (10% test, 90% train)
set.seed(42) # Setting seed for reproducibility
cells_90 <- sample(all_cells, size = length(all_cells) * 0.9)

# Split the WT object into two
WT_90 <- subset(WT, cells = cells_90)
WT_10 <- subset(WT, cells = setdiff(all_cells, cells_90))

# PCA
#####################################################################################################################################

# TRAIN:
# Reference method using PCA, binary (present or absent) and median cell state proportions:

train <- WT_90@meta.data
train <- train %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))

all_combinations <- expand.grid(embryo_stage = unique(train$embryo_stage), 
                                cluster_name = cluster_names)


# calculate median cell state poportions 

total_cells <- train %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
proportions <- train %>%
  left_join(total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

# To get PC1 distribution of cell state proportions in each train embryo 
proportions_complete <- all_combinations %>%
  left_join(proportions, by = c("cluster_name" = "cluster_names", "embryo_stage")) %>%
  replace_na(list(count = 0, proportion = 0))

# Calculate median proportions for every stage and cluster_name over embryos
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

# Prepare metadata for join
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

# Done: Binary and Median proportions 

# PCA for each stage, on each embryo (median_proportion vector and on binary info vectors)
# Reshape and perform PCA on median proportions
median_proportions_wide <- median_proportions %>% select(-median_proportion) %>%
  pivot_wider(names_from = cluster_name, values_from = normalized_median_proportion, values_fill = list(normalized_median_proportion = 0.0)) 
  
# Identify columns with near zero variance
nzv_prop <- nearZeroVar(median_proportions_wide)
# "38_PostProxPrimStreak" 

# Remove near zero variance columns from the dataset (if there are nzvs)
median_proportions_wide_clean <- median_proportions_wide %>% select(-nzv_prop)
median_proportions_wide_clean <- median_proportions_wide_clean %>% select(-stage)

# PCA on median proportions
pca_median_proportions <- prcomp(median_proportions_wide_clean, retx=TRUE, center = TRUE, scale. = TRUE)

# Convert PCA results to a data frame
pca_prop_df <- as.data.frame(pca_median_proportions$x)

# Add the stage name as a new column
pca_prop_df$stage <- stages

# Reshape and perform PCA on binary presence
binary_presence_wide <- median_binary %>%
  pivot_wider(names_from = cluster_name, values_from = median_binary, values_fill = list(median_binary = 0)) %>%
  select(-stage) 

# Remove near zero variance columns from the dataset
nzv_bin <- nearZeroVar(binary_presence_wide)
#  9 22 32 34 39
binary_presence_wide_clean <- binary_presence_wide[, -nzv_bin]

pca_binary_presence <- prcomp(binary_presence_wide_clean, retx=TRUE, center = TRUE, scale. = TRUE)

# Convert them into dataframes
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

# Calculate proportion of variance explained by each component 
proportion_variance_bin <- (variance_bin / total_variance_bin)*100
proportion_variance_prop <- (variance_prop / total_variance_prop)*100

unique_stages <- unique(WT@meta.data$stage)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, unique_stages)
combined_df$e_days <- e_days

# colors for each stage
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
ggsave("autoplot_pca_binary_presence.pdf", autoplot(pca_binary_presence))

# Scatter-Plots for 2 principal components 
scatter_bin <- ggplot(df_binary_presence, aes(x=PC1_Binary, y=PC2_Binary, label=e_days))+
  geom_point()+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC2 binary (", round(proportion_variance_bin[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  geom_text(hjust=0.5, vjust=2)
ggsave("PC1_PC2_values_median_WT_binary.pdf", scatter_bin)

scatter_prop <- ggplot(df_median_proportions, aes(x=PC1_Proportion, y=PC2_Proportion, label=e_days))+
  geom_point()+
  xlab(paste0("PC1 proportions (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC2 proportions (", round(proportion_variance_prop[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  geom_text(hjust=0.5, vjust=2)
ggsave("PC1_PC2_values_normalized_median_WT_proportions.pdf", scatter_prop)

#####################################################################################################################

# Distributions of principal components for each embryo 

binary_presence$embryo_stage <- paste(binary_presence$embryo,binary_presence$stage, sep="__")

# Reshape and perform PCA on binary presence
binary_presence_wide2 <- binary_presence %>%
  group_by(embryo_stage) %>%
  pivot_wider(names_from = cluster_name, values_from = binary_presence, values_fill = list(binary_presence = 0))

proportions_complete <- proportions_complete %>%
  mutate(
    embryo = stringr::str_extract(embryo_stage, "embryo\\d+"), # Extract 'embryo' followed by digits
    stage = stringr::str_extract(embryo_stage, "WT_\\d+") 
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


# PCA on all train embryos 

# prediction of PCs for embryos of train WT_90 -> loadings / rotations of WT_90 (train medians) were used 
pca_train_proportions <- predict(pca_median_proportions, newdata=proportions_wide2_clean)
pca_train_binary <-  predict(pca_binary_presence, newdata=binary_wide2_clean)  

# Convert PCA results to a data frame
pca_prop_df2 <- as.data.frame(pca_train_proportions)
pca_bin_df2 <- as.data.frame(pca_train_binary)


# Add embryo and stage information back to it 
pca_prop_df2$embryo_stage <- embryo_stage_prop
pca_prop_df2$stage <- stages_prop
pca_bin_df2$embryo_stage <- embryo_stage_bin

# Convert them into dataframes
train_combined_df <- 
  merge(data.frame(embryo_stage=pca_prop_df2$embryo_stage, PC1_Proportion = pca_prop_df2[,"PC1"], stage =pca_prop_df2$stage),data.frame(embryo_stage = pca_bin_df2$embryo_stage, PC1_Binary = pca_bin_df2[,"PC1"]), by="embryo_stage")%>%
  mutate(e_day = named_stages[stage]) 

# PC1 Distribution Plots: 
# Plot for PC1_Proportion
pc1_prop <- ggplot(train_combined_df, aes(x = e_day, y = PC1_Proportion)) +
  geom_boxplot() +
  theme_classic() +
  labs( x = "Stage", y = "PC1 Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for PC1_Binary
pc1_bin <- ggplot(train_combined_df, aes(x = e_day, y = PC1_Binary)) +
  geom_boxplot() +
  theme_classic() + 
  labs(x = "Stage", y = "PC1 Binary") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("PC1_distr_prop.pdf",width=8, pc1_prop)
ggsave("PC1_distr_bin.pdf",width=8, pc1_bin)

################################################################################################################################################################################

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

# stacked plots with WT as reference 
median_stacked_plot <- ggplot(median_proportions, aes(x = e_day, y = normalized_median_proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Stage", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45) )
ggsave("/project/PRCcomp/images/Developmental_Delay/median_stacked_cellstates_normalized.pdf",width=8, median_stacked_plot)

median_binary$combined_cluster <- paste(median_binary$cluster_name, median_binary$cluster_germlayer, sep = "_")
median_binary$combined_cluster <- factor(median_binary$combined_cluster, levels = unique(median_binary$combined_cluster))

# binary tile plots 
binary_plot <- ggplot(median_binary, aes(x = e_day, y = cluster_name, fill = factor(median_binary))) +
  geom_tile(color = "grey") + # Add tiles with a grey border
  scale_fill_manual(values = c("1" = "black", "0" = "white", "0.5" = "grey"),  labels = c("0" = "0 (absent)", "1" = "1 (present)")) + # Set fill colors
  theme_minimal() + # Use a minimal theme
  facet_grid(cluster_germlayer ~ ., scales = "free_y", space = "free_y")+ # Group by germlayer with space 
  labs(x = "Stage", y = "Cell State", fill = "Median Binary") + # Add labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) # Remove minor grid lines
ggsave("median_binary_plot.pdf",width=8, height=14, binary_plot)


##############################################################################################################################################################
# VALIDATION:
##############################################################################################################################################################

test <- WT_10@meta.data
test <- test %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))

# same amount of embryos (72) as in test (all are included):
length(unique(test$embryo_stage))
length(unique(train$embryo_stage))

# and similar stage proportions: 
stage_proportions_test <- test %>% group_by(stage) %>% summarise(proportion = n()/dim(test)[1])
stage_proportions_train <- train %>% group_by(stage) %>% summarise(proportion = n()/dim(train)[1])

test_all_combinations <- expand.grid(embryo_stage = unique(test$embryo_stage), 
                                cluster_name = cluster_names)

# calculate median cell state poportions 
test_total_cells <- test %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
test_proportions <- test %>%
  left_join(test_total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

test_proportions_complete <- test_all_combinations %>%
  left_join(test_proportions, by = c("cluster_name" = "cluster_names", "embryo_stage")) %>%
  replace_na(list(count = 0, proportion = 0)) %>%
  mutate(split=embryo_stage)%>%
  separate(split, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo_stage, embryo, stage, cluster_name, proportion = proportion)

test_proportions_wide <- test_proportions_complete %>%
  pivot_wider(names_from = cluster_name, values_from = proportion, values_fill = list(proportion = 0.0)) 

test_embryo_stage <- test_proportions_wide$embryo_stage

# Extract column names from median_proportions_wide
column_order <- colnames(median_proportions_wide)

# Reorder columns of test_proportions_wide to match
test_proportions_wide <- test_proportions_wide %>% 
  select(all_of(column_order)) 
 
test_proportions_wide_clean <- test_proportions_wide %>% select(-stage)

# prediction of PCs for test dataset WT_10 -> loadings / rotations of WT_90 were used 
pred_prop <- predict(pca_median_proportions, newdata=test_proportions_wide_clean)  
# Convert PCA results to a data frame
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

# stacked plots with WT as reference 
test_stacked_plot <- ggplot(test_proportions_complete, aes(x = factor(order_number), y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        facet_grid(~e_day, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45) )
ggsave("test_stacked_cellstates.pdf",width=14, test_stacked_plot)


# Binary Information of cell states scatter

# Prepare metadata for join
test_prepared <- test %>%
  distinct(embryo_stage, cluster_names) %>%
  mutate(presence = 1)

# Merge with metadata to check presence
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
        mutate(order_number = (match(embryo_stage, unique(embryo_stage))))

# binary tile plots 
test_binary_plot <- ggplot(test_binary_presence, aes(x = factor(order_number), y = cluster_name, fill = factor(binary_presence))) +
  geom_tile(color = "grey") + # Add tiles with a grey border
  scale_fill_manual(values = c("1" = "black", "0" = "white"),  labels = c("0" = "0 (absent)", "1" = "1 (present)")) + # Set fill colors
  theme_minimal() + # Use a minimal theme
  facet_grid(cluster_germlayer ~  e_day , scales = "free", space = "free_y")+ # Group by germlayer with space 
  labs(x = "Embryo", y = "Cell State", fill = "") + # Add labels
  theme(axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels for better readability
        panel.grid.major = element_blank(), # Remove major grid lines
        panel.grid.minor = element_blank()) # Remove minor grid lines
ggsave("test_binary_plot_b.pdf",width=8, height=14, test_binary_plot)

test_binary_wide <- test_binary_presence %>%
  select(-c(order_number, e_day, cluster_germlayer, stage, embryo)) %>%
  pivot_wider(names_from = cluster_name, values_from = binary_presence) 

test_embryo_stage_bin <- test_binary_wide$embryo_stage
test_binary_wide <- select(test_binary_wide,-embryo_stage)

# Extract column names from median_proportions_wide
column_order <- colnames(binary_presence_wide_clean)

# Reorder columns of test_proportions_wide to match
test_binary_wide <- select(test_binary_wide,all_of(column_order)) 

# prediction of PCs for test dataset WT_10 -> loadings / rotations of WT_90 were used 
pred_bin <- predict(pca_binary_presence, newdata=test_binary_wide)  
# Convert PCA results to a data frame
pred_bin <- as.data.frame(pred_bin)

# Add embryo and stage information back to it 
pred_bin$embryo_stage <- test_embryo_stage_bin

# Convert them into dataframes
test_combined_df <- 
merge(data.frame(embryo_stage=pred_prop$embryo_stage, PC1_Proportion = pred_prop[,"PC1"]),data.frame(embryo_stage = pred_bin$embryo_stage, PC1_Binary = pred_bin[,"PC1"]), by="embryo_stage")

test_combined_df <- test_combined_df %>%
  mutate(split=embryo_stage)%>%
  separate(split, into = c("embryo", "stage"), sep = "__") %>%
  select(embryo_stage, embryo, stage, PC1_Proportion,PC1_Binary)

# Add a column for named stages
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
ggsave("PC1_values_median_WT_embryos_a.pdf",width=8, scatter1)

common_cols <- intersect(colnames(combined_df), colnames(test_combined_df))

# Subset both dataframes to only include common columns
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
ggsave("PC1_values_median_WT_embryos_b.pdf",width=8, scatter2)

test_scatter <- ggplot(test_combined_df, aes(x=PC1_Binary, y=PC1_Proportion, color=e_days))+
  geom_point(aes(colour=e_days,shape=Group))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  xlab(paste0("PC1 binary (", round(proportion_variance_bin[1]), "%)"))+
  ylab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  theme_minimal()+
  theme_classic()+
  labs(colour="Stage")
ggsave("PC1_values_test_WT_embryos_a.pdf", test_scatter+
  geom_text(aes(label=e_days),hjust=0.5, vjust=2))

together <- test_scatter+ geom_point(data=combined_df, aes(x=PC1_Binary, y=PC1_Proportion, , shape=Group)) + 
  geom_text(data=combined_df, aes(label=e_days),hjust=0.5, vjust=2)
ggsave("PC1_values_test_median_WT_embryos_b.pdf",together)


# For marginal boxplots:

train_combined_df <- train_combined_df %>% rename("e_days"="e_day")
common_cols <- intersect(colnames(combined_full_df), colnames(train_combined_df))
# Combine the dataframes
combined_full_df <- rbind(combined_full_df[, common_cols] , train_combined_df[, common_cols])

scatter3 <- scatter2 + 
  geom_point(data = subset(combined_full_df, Group == "train (median)"), 
             aes(shape=Group, color=Group), size=3, color="black") +
  geom_text(data = subset(combined_full_df, Group == "train (median)"),
            aes(label=e_days), hjust=0.5, vjust=2, check_overlap = TRUE, color="black")+
  ylim(-10.5,7.5)

l2 <- get_legend(scatter3)
scatter3 <- scatter3 + theme(legend.position="none") 
ggsave("PC1_values_test_median_WT_embryos_black.pdf",scatter3)
subset <- subset(combined_full_df, Group == "train" | Group == "test")            
# Marginal boxplot of x (top panel) and y (right panel)
#xplot <- ggMarginal(scatter3,data = subset, x = "e_days", y = "PC1_Binary",type="boxplot", color=Group, fill=Group, margins="x")

# Marginal boxplot of x (top panel) and y (right panel)
xplot <- ggboxplot(subset, x = "e_days", y = "PC1_Binary", 
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  xlab("Stage")+
  rotate()+
  rremove("xlab")+
  ylim(-6.5,6.5)+
  theme(axis.text.x = element_blank())
yplot <- ggboxplot(subset, x = "e_days", y = "PC1_Proportion",
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  rremove("legend")+
  rremove("ylab")+
  xlab("Stage")+
  ylim(-10.5,7.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank())
l1 <- get_legend(xplot)
xplot <- xplot+rremove("legend")
complete <- plot_grid(xplot, NULL, scatter3, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2, 3))
ggsave("PC1_values_scatter_boxplots.pdf",complete)
ggsave("PC1_values_scatter_boxplots_legend.pdf",ggarrange(l1,l2))

p1 <- insert_xaxis_grob(scatter3, xplot, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, yplot, grid::unit(.2, "null"), position = "right")
ggsave("PC1_values_scatter_boxplots_b.pdf",p2)

####################################################################################################################################################################
# Done: PCA on test
# Now: classification by euclidean distance of PC1s to train PC1s 

euclidean_distance <- function(point1, point2) {
  if (length(point1) != length(point2)) {
    stop("Points must have the same dimensions")
  }
  return(sqrt(sum((point1 - point2)^2)))
}
# Initialize columns in test_combined_df for each stage
for (stage in unique_stages) {
  test_combined_df[[stage]] <- NA
}
# Loop through each row in test_combined_df
for (i in 1:nrow(test_combined_df)) {
  point_test <- test_combined_df[i, c("PC1_Binary", "PC1_Proportion")]
  # Compute distances to each point in combined_df
  for (j in 1:nrow(combined_df)) {
    point_combined <- combined_df[j, c("PC1_Binary", "PC1_Proportion")]
    distance <- euclidean_distance(point_test, point_combined)
    stage_name <- combined_df$Stage[j]
    test_combined_df[i, stage_name] <- distance
  }
}
test_combined_df$closest_stage <- NA
# Loop through each row of test_combined_df
for (i in 1:nrow(test_combined_df)) {
  # Extract the distances for the current row, excluding non-distance columns
  distances <- test_combined_df[i, unique_stages]
  # Find the minimum distance and the corresponding stage name
  min_distance <- min(distances, na.rm = TRUE)
  closest_stage <- names(distances)[which.min(distances)]
  # Update the result column with the name of the closest stage
  test_combined_df$predicted_stage[i] <- closest_stage
}

# Classification rate for stages
correct_predictions <- sum(test_combined_df$stage == test_combined_df$predicted_stage)
total_predictions <- nrow(test_combined_df)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (PC 1 of binary pressence and proportions of cell states):", classification_rate))


###############################################################################################################################
# Classification on PC1 and PC2 of cell state proportions (wihout binary)
###############################################################################################################################

# df_medain_proportions: PC1 and PC2 values of median (train)
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

common_cols <- intersect(colnames(df_train_pca), colnames(df_test_pca))
# Combine the dataframes
df_combi_pca <- rbind( df_train_pca[, common_cols], df_test_pca[, common_cols])

scatter_prop <- scatter_prop +  coord_cartesian(ylim = c(-10.5, 5), xlim = c(-9, 9))

# Ensure 'e_day' is a factor and reverse its levels
#df_combi_pca$e_day <- factor(df_combi_pca$e_day, levels = rev(unique(df_combi_pca$e_day)))

# Marginal boxplot of x (top panel) and y (right panel)
xplot <- ggboxplot(df_combi_pca, x = "e_day", y = "PC1_Proportion", 
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  xlab("Stage")+
  rotate()+
  rremove("xlab")+
  ylim(-9,9)+
  theme(axis.text.x = element_blank())+
  scale_x_discrete(limits=rev)
yplot <- ggboxplot(df_combi_pca, x = "e_day", y = "PC2_Proportion",
  color = "Group", fill = "Group", ggtheme = theme_classic(), alpha=0.5)+
  rremove("legend")+
  rremove("ylab")+
  xlab("Stage")+
  ylim(-10.5,5)+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_blank())+
  scale_x_discrete(limits=rev)
l1 <- get_legend(xplot)
xplot <- xplot+rremove("legend")
complete <- plot_grid(xplot, NULL, scatter_prop, yplot, ncol = 2, align = "hv", 
          rel_widths = c(3, 2), rel_heights = c(2, 3))
ggsave("PC1_PC2_proportion_values_scatter_boxplots.pdf",complete)
ggsave("PC1_PC2_proportion_values_scatter_boxplots_legend.pdf",l1)

test_scatter <- ggplot(df_test_pca, aes(x=PC1_Proportion, y=PC2_Proportion, color=e_day))+
  geom_point(aes(colour=e_day,shape=Group))+
  scale_shape_manual(values=c("train (median)"=8, "test"=16))+
  scale_colour_brewer(palette="Dark2")+
  xlab(paste0("PC1 proportions (", round(proportion_variance_prop[1]), "%)"))+
  ylab(paste0("PC2 proportions (", round(proportion_variance_prop[2]), "%)"))+
  theme_minimal()+
  theme_classic()+
  labs(colour="Stage")
ggsave("PC1_2_values_test_WT_embryos.pdf", test_scatter)

colnames(df_median_proportions)[colnames(df_median_proportions) == "e_days"] <- "e_day"
together <- test_scatter + 
  geom_point(data=df_median_proportions, aes(x=PC1_Proportion, y=PC2_Proportion, shape=Group)) + 
  geom_text(data=df_median_proportions, aes(label=e_day), hjust=0.5, vjust=2)
ggsave("PC1_2_test_median_WT_embryos.pdf",together)


# classification by euclidean distance of PC1 & PC2 to median (train) PC1 & PC2  

# Initialize columns in test for each stage
for (stage in e_days) {
  df_test_pca[[stage]] <- NA
}
# Loop through each row in test
for (i in 1:nrow(df_test_pca)) {
  # Get the current point from test
  point_test <- df_test_pca[i, c("PC1_Proportion", "PC2_Proportion")]

  # Compute distances to each point in train 
  for (j in 1:nrow(df_median_proportions)) {
    point_train <- df_median_proportions[j, c("PC1_Proportion", "PC2_Proportion")]
    distance <- euclidean_distance(point_test, point_train)
    #print(distance)
    stage_name <- df_median_proportions$e_day[j]
    df_test_pca[i, stage_name] <- distance
  }
}
df_test_pca$closest_stage <- NA

# Loop through each row of test df 
for (i in 1:nrow(df_test_pca)) {
  # Extract the distances for the current row, excluding non-distance columns
  distances <- df_test_pca[i, e_days]
  # Find minimum distance and the corresponding stage name
  min_distance <- min(distances, na.rm = TRUE)
  closest_stage <- names(distances)[which.min(distances)]
  #print(closest_stage)
  # Update the result column with the name of the closest stage
  df_test_pca$predicted_stage[i] <- closest_stage
}

# Classification rate for stages
correct_predictions <- sum(df_test_pca$e_day == df_test_pca$predicted_stage)
total_predictions <- nrow(df_test_pca)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (PC 1 & 2 of cell states proportions):", classification_rate))

# Classification Statistics (confusion matrix included)
conf_mat <- confusionMatrix(table(df_test_pca$predicted_stage, df_test_pca$e_day))
print(conf_mat)

# Confusion Matrix:
confMat_melt <- as.data.frame(as.table(conf_mat$table/colSums(conf_mat$table)))
confMat_plot <- ggplot(data = confMat_melt, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Prediction", y = "Stage", fill = "Normalized Frequency")
ggsave("PC1_PC2_proportion_confusion.pdf",width=10, confMat_plot)

actual_numeric <- as.double(substr(df_test_pca$e_day, 2, 4))
predicted_numeric <- as.double(substr(df_test_pca$predicted_stage, 2, 4))

mae <- mae(actual_numeric, predicted_numeric)
mse <- mse(actual_numeric, predicted_numeric)

print(paste0("RMSE = ", sqrt(mse)))

# test_proportions_complete: all cell state proportions for test embryos
# df_test_pca: all PC1 & 2 values for all test embryo cells and their real and predicted stage (e_day)
df_test_pca <- df_test_pca %>% mutate(PC_angle = 180*atan2(PC1_Proportion, PC2_Proportion)/pi)
df_test_pca <- df_test_pca %>% arrange(predicted_stage,desc(PC_angle))
df_test_pca <- df_test_pca %>% mutate(embryo_stage = paste(embryo, stage, sep = "__"))

test_proportions_complete <- left_join(test_proportions_complete, 
                                       select(df_test_pca, embryo_stage,predicted_stage, PC_angle), 
                                       by = "embryo_stage", 
                                       copy = TRUE) # Add copy = TRUE to address the src issue

test_proportions_complete <- test_proportions_complete %>%
  arrange(predicted_stage, desc(PC_angle)) 

embryo_stages_order <- test_proportions_complete %>%
  distinct(embryo_stage, .keep_all = TRUE) %>%
  mutate(order_number = row_number()) %>%
  select(embryo_stage, order_number)

test_proportions_complete <- test_proportions_complete %>%
  select(-order_number)%>%
  left_join(embryo_stages_order, by = "embryo_stage")

#### To order embryo by assigned stage (angle of (PC1,PC2) vector)

# Stacked Plot for each test embryo: 
test_proportions_complete <- test_proportions_complete %>%
  arrange(cluster_germlayer, cluster_name)
test_proportions_complete$cluster_name <- factor(test_proportions_complete$cluster_name, levels = unique(test_proportions_complete$cluster_name))
test_proportions_complete$x_labels <- paste0(test_proportions_complete$order_number,"_",round(test_proportions_complete$PC_angle, 3),"_", test_proportions_complete$e_day)
# Assuming test_proportions_complete$order_number is numeric or character
test_proportions_complete$order_number <- factor(test_proportions_complete$order_number, levels = unique(test_proportions_complete$order_number))

# stacked plots with WT as reference 
test_stacked_plot <- ggplot(test_proportions_complete, aes(x = factor(x_labels), y = proportion*100, fill = cluster_name)) +
        geom_bar(stat = "identity", position = "stack") +
        facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        facet_grid(~predicted_stage, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45) )
ggsave("test_stacked_cellstates_prediction_order.pdf",width=14, test_stacked_plot)

# For the comparison of PC1 bin+prop classification 
# test_combined_df: PC1 Binary and Proportion and actual and predicted stage 
test_combined_df <- test_combined_df %>% mutate(predicted_e_day = named_stages[predicted_stage])

# Classification Statistics (confusion matrix included)
conf_mat <- confusionMatrix(table(test_combined_df$predicted_e_day, test_combined_df$e_days))
print(conf_mat)

# Confusion Matrix:
confMat_melt <- as.data.frame(as.table(conf_mat$table/colSums(conf_mat$table)))
confMat_plot <- ggplot(data = confMat_melt, aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Prediction", y = "Stage", fill = "Normalized Frequency")
ggsave("PC1_PC1_prop_bin_confusion.pdf",width=10, confMat_plot)

actual_numeric <- as.double(substr(test_combined_df$e_days, 2, 4))
predicted_numeric <- as.double(substr(test_combined_df$predicted_e_day, 2, 4))

mae <- mae(actual_numeric, predicted_numeric)
mse <- mse(actual_numeric, predicted_numeric)

print(paste0("RMSE = ", sqrt(mse)))

###############################################################################################################################################################

# training set:
save(WT_90, file="WT_90.Robj")
# validation set:
save(WT_10, file="WT_10.Robj")
# PCA models: 
save(pca_binary_presence, file="pca_binary_presence.Robj")
save(pca_median_proportions, file="pca_median_proportions.Robj")



