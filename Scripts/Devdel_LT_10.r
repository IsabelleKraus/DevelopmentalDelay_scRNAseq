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
# Define a set of distinct colors for each stage (same as in PCA)
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames(stage_colors, e_days)
color_palette <- setNames(stage_colors, c("6.5","7","7.5","8","8.5","9"))

# Determine the size of subsets (10% test, 90% train)
set.seed(42) # Setting seed for reproducibility
cells_90 <- sample(all_cells, size = length(all_cells) * 0.9)

# training:
WT_90 <- subset(WT, cells = cells_90)
# validation:
WT_10 <- subset(WT, cells = setdiff(all_cells, cells_90))

red_palette <- c(
  "6.5" = "#FFCCCC",  # Lighter red
  "7" = "#FF9999",
  "7.5" = "#FF6666",
  "8" = "#FF3333",
  "8.5" = "#FF0000",
  "9" = "#CC0000"   # Darker red
)

# Label Transfer:
#################################################################################################
# Parameters:
red <- "pcaproject"
kfilter <- NA

# WT_90 as training reference, WT_10 as test query
plan("multicore", workers = 2)

# transfer stages 
anchors <- FindTransferAnchors(reference = WT_90, query = WT_10, k.filter = kfilter, reduction = red)    
prediction1 <- TransferData(anchorset = anchors, refdata = WT_90@meta.data$stage, weight.reduction = red)

WT_10 <- AddMetaData(object = WT_10, metadata = prediction1)

pred <- WT_10@meta.data$predicted.id
real_stage <- WT_10@meta.data$stage
rownames_df <- rownames(WT_10@meta.data)

df <- data.frame(stage = real_stage, predicted_stage = pred)
rownames(df) <- rownames_df

prediction1$real_stage <- real_stage

# Calculate mean prediction scores for each real_stage
mean_scores_by_stage <- prediction1 %>%
  group_by(real_stage) %>%
  summarise(across(starts_with("prediction.score"), mean, na.rm = TRUE))

mean_scores_by_stage <- mean_scores_by_stage %>%
  select(-last_col())

mean_scores_by_stage$real_stage <- e_days  

# Prepare the data for plotting
mean_scores_long <- mean_scores_by_stage %>%
  pivot_longer(-real_stage, names_to = "score_type", values_to = "score")

pred_e <- c("prediction.score.WT_90"="E9.0","prediction.score.WT_85"="E8.5","prediction.score.WT_80"="E8.0","prediction.score.WT_75"="E7.5","prediction.score.WT_70"="E7.0","prediction.score.WT_65"="E6.5")

hmap <- ggplot(mean_scores_long, aes(y = real_stage, x = score_type, fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels=pred_e)+
  theme_minimal() +
  labs(y = "Stage", x = "Prediction", fill = "Mean score") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("hmap_LT_pred_10.pdf",width=8, hmap)

WT_10$e_day <- named_stages[WT_10$stage]

# classification: stage with highest prediction score
WT_10$pred_e <- named_stages[WT_10$predicted.id]

# UMAPs:
umap_pred <- DimPlot(WT_10, cols=color_mapping, group.by="pred_e")
ggsave("umap_LT_pred_10.pdf",width=8, umap_pred)

umap_real <- DimPlot(WT_10, cols=color_mapping, group.by="e_day")
ggsave("umap_LT_actual_stages_10.pdf",width=8, umap_real)

correct <- WT_10[,WT_10@meta.data$stage == WT_10@meta.data$predicted.id]
false <- WT_10[,WT_10@meta.data$stage != WT_10@meta.data$predicted.id]
umap_correct <- DimPlot(correct, cols=color_mapping, group.by="e_day")
umap_false <- DimPlot(false, cols=color_mapping, group.by="pred_e")

ggsave("umap_LT_correct_predictions_10.pdf",width=8,umap_correct)
ggsave("umap_LT_false_predictions_10.pdf",width=8, umap_false)

# deltas of false predictions:
false$predicted_numeric <- as.double(paste0(substr(false$predicted.id, 4, 4),".",substr(false$predicted.id, 5, 5)))
false$actual_numeric  <- as.double(paste0(substr(false$stage, 4, 4),".",substr(false$stage, 5, 5)))
false$delta_pred <- false$predicted_numeric - false$actual_numeric

myColorPalette <- colorRampPalette(c("blue","grey", "red"))
colors <- myColorPalette(11)
# colors for 0.5 deltas
delta_colors <- setNames(colors, seq(-2.5,2.5,0.5))

umap_delta <- DimPlot(false,cols=delta_colors,  group.by="delta_pred")
ggsave("umap_LT_delta_false_predictions_10.pdf",width=8,umap_delta)

all_WT_10 <- WT_10
all_WT_10$predicted_numeric <- as.double(paste0(substr(all_WT_10$predicted.id, 4, 4),".",substr(all_WT_10$predicted.id, 5, 5)))
all_WT_10$actual_numeric  <- as.double(paste0(substr(all_WT_10$stage, 4, 4),".",substr(all_WT_10$stage, 5, 5)))
all_WT_10$delta_pred <- all_WT_10$predicted_numeric - all_WT_10$actual_numeric

umap_delta_all <- DimPlot(all_WT_10,cols=delta_colors,  group.by="delta_pred")
ggsave("umap_LT_delta_predictions_10.pdf",width=8,umap_delta_all)

# grey WT_90 below those 
umap_data <- as.data.frame(WT_90[["umap"]]@cell.embeddings)
umap_data$State <- WT@meta.data[row.names(umap_data),"State"]
umap_data$stage <- WT@meta.data[rownames(umap_data), "stage"]
umap_data <- umap_data %>% left_join(annotation, by=c("State"="State"))
head(umap_data)

# UMAP Plot per Stage with all stages underlying in grey:
umap_data_stage <- dplyr::select(umap_data, -stage) 
stage.labs = e_days
names(stage.labs) <- unique_stages
umap_data$e_day <- named_stages[umap_data$stage]

facet_plot <- ggplot(umap_data, aes(x= UMAP_1, y=UMAP_2)) +
  geom_point(data=umap_data_stage, size=0.5, colour="grey82", alpha=0.5) +
  theme_minimal() +
  scale_color_manual(name="Stage", values = color_mapping) +
  facet_wrap(~ e_day, ncol = length(unique_stages), strip.position="top")+
  geom_point(aes(color=e_day), size=0.5) + 
  theme_void() +
  theme(strip.text=element_text(size=20), panel.border = element_rect(colour = "black", linewidth=1))+
  guides(color="none")
ggsave("WT_UMAP_stages_vs_all_wrapped_plot_LT_10.pdf",facet_plot, width=6*length(unique_stages), height=5)

# Classification rate 
correct_predictions <- sum(df$stage == df$predicted_stage)
total_predictions <- nrow(df)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages (red=",red,", kfilter=",kfilter,"): ", classification_rate))

###########################################################################################
# Classification for whole embryos using (all) prediction scores per cell 
# Per Embryo: Mean stage of weighted means (cell stages)

predicted_numeric <- as.double(paste0(substr(WT_10@meta.data$predicted.id, 4, 4),".",substr(WT_10@meta.data$predicted.id, 5, 5)))
actual_numeric <- as.double(paste0(substr(WT_10@meta.data$stage, 4, 4),".",substr(WT_10@meta.data$stage, 5, 5) ))

# test data:
WT_10 <- AddMetaData(object = WT_10, metadata = predicted_numeric, col.name="predicted_numeric")
WT_10 <- AddMetaData(object = WT_10, metadata = actual_numeric, col.name="actual_numeric")
data <- WT_10@meta.data

prediction_scores <- data[, grep("prediction.score.WT_", names(data))]
stages_numeric <- c(9.0,8.5,7.0,8.0,7.5,6.5)
# Compute weighted sum of stages
weighted_sums <- apply(prediction_scores, 1, function(x) sum(x * stages_numeric))
data$weighted_sums <- weighted_sums

# Group by stage and calculate for each embryo within stage the predicted_numeric stage:
means <- data %>% group_by(actual_numeric, embryo) %>%  summarise(mean = mean(weighted_sums)) %>%
    ungroup()
means$closest_stage <- sapply(means$mean, function(x) {
  stages_numeric[which.min(abs(x - stages_numeric))]
})

# classified test data 
# Merge the 'means' dataframe back into the 'data' dataframe based on 'actual_numeric' and 'embryo'
merged_data <- merge(data, means, by = c("actual_numeric", "embryo"), all.x = TRUE)
merged_data <- merged_data %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))

test_all_combinations <- expand.grid(embryo_stage = unique(merged_data$embryo_stage), 
                                cluster_name = cluster_names)
                        
test_total_cells <- merged_data %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
test_proportions <- merged_data %>%
  left_join(test_total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total ) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

merged_data_updated <- merged_data %>%
  left_join(test_proportions %>% select(cluster_names, embryo_stage, proportion), 
            by = c("cluster_names", "embryo_stage")) %>%
  mutate(proportion = if_else(is.na(proportion), 0, proportion))
head(merged_data_updated)

grouped_df <- merged_data_updated %>%
  group_by(embryo_stage, cluster_names) %>%
  summarise(
    proportion = first(proportion), 
    mean = first(mean),
    closest_stage = first(closest_stage),
    actual_numeric = first(actual_numeric),
    cluster_germlayer = first(cluster_germlayer)
  ) %>%
  ungroup()  
head(grouped_df)

sum <- grouped_df %>% group_by(embryo_stage) %>% summarize(wm=sum(proportion))

grouped_df <- grouped_df %>%
  arrange(cluster_germlayer, cluster_names) 
grouped_df$cluster_names <- factor(grouped_df$cluster_names, levels = unique(grouped_df$cluster_names))

# stacked plots with WT as reference 
test_stacked_plot <- ggplot(grouped_df, aes(x = factor(round(mean,5)), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo Stage Prediction", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~closest_stage, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))

ggsave("test_stacked_cellstates_ordered_by_prediction_stage_means_LT_10.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
grouped_df$delta <- grouped_df$mean - grouped_df$actual_numeric
delta_df <- grouped_df %>% group_by(embryo_stage,actual_numeric,closest_stage, mean, delta) %>%  summarize()

p <- ggplot(delta_df, aes(x = factor(round(mean,5)), fill = delta)) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0,limits = c(-2.5, 2.5)) + # Color transition
  labs(x = "", y = "", fill = "Predction-Actual") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lp <- get_legend(p)
ggsave("prediction_delta_LT_10.pdf",width=21,height=2, p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo Stage Prediction",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())

# To align delta and stacked plot with actual (true) heatmap 
a <- ggplot(delta_df, aes(x = factor(round(mean,5)), fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual Stage", values=color_palette) +
  labs(x = "", y = "", fill = "Actual Stage") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
la <- get_legend(a)

delta_df$classification <- ifelse(delta_df$actual_numeric==delta_df$closest_stage,"true","false")

b <- ggplot(delta_df, aes(x = factor(round(mean,5)), fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("false"="red","true"="olivedrab")) +
  labs(x = "", y = "", fill = "Classification") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 

lb <- get_legend(b)

a <- a + rremove("xlab") + rremove("legend")
p <- p + rremove("xlab") + rremove("legend")
b <- b + rremove("xlab") + rremove("legend")
complete <- b/a/p/test_stacked_plot + plot_layout(heights = c(0.2,0.2,0.2,4))
ggsave(paste0("test_stacked_delta_actual_classification_complete_LT_10.pdf"),width=23,height=9,complete)
ggsave(paste0("test_stacked_delta_actual_classification_complete_legend_LT_10.pdf"),ggarrange(lb,lp,la))


# Confusion Matrix:
conf_mat <- confusionMatrix(factor(means$closest_stage), factor(means$actual_numeric))
print(conf_mat)

confMat_melt <- as.data.frame(as.table(conf_mat$table)/colSums(conf_mat$table))
confMat_melt$Prediction <- factor(confMat_melt$Prediction)
confMat_melt$Reference <- factor(confMat_melt$Reference)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]

named_stages2 <- c("6.5" = "E6.5","7" = "E7.0", "7.5" = "E7.5","8" = "E8.0", "8.5" = "E8.5","9" = "E9.0")

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Prediction), y =factor(Reference), fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels=named_stages2)+
  scale_y_discrete(labels=named_stages2)+
  theme_minimal() +
  labs(x = "Prediction", y = "Stage", fill = "Normalized Frequency")
ggsave("LT_10_all_pred_scores_confusion.pdf",width=10, confMat_plot)

mae <- mae(means$actual_numeric, means$mean)
mse <- mse(means$actual_numeric, means$mean)
print(paste0("RMSE = ", sqrt(mse)))

mae <- mae(means$actual_numeric, means$closest_stage)
mse <- mse(means$actual_numeric, means$closest_stage)
print(paste0("RMSE = ", sqrt(mse)))

##############################################################################################################################################

# Group by stage and calculate for each embryo within stage the weighted  predicted_numeric stage:
##############################################################################################################################################

means <- data %>% group_by(actual_numeric, embryo) %>%  summarise(weighted_mean = weighted.mean(predicted_numeric, prediction.score.max)) %>%
    ungroup()
means$closest_stage <- sapply(means$weighted_mean, function(x) {
  stages_numeric[which.min(abs(x - stages_numeric))]
})

# classified test data 
# Merge the 'means' dataframe back into the 'data' dataframe based on 'actual_numeric' and 'embryo'
merged_data <- merge(data, means, by = c("actual_numeric", "embryo"), all.x = TRUE)
merged_data <- merged_data %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))

test_all_combinations <- expand.grid(embryo_stage = unique(merged_data$embryo_stage), 
                                cluster_name = cluster_names)                  

test_total_cells <- merged_data %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
test_proportions <- merged_data %>%
  left_join(test_total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total ) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

merged_data_updated <- merged_data %>%
  left_join(test_proportions %>% select(cluster_names, embryo_stage, proportion), 
            by = c("cluster_names", "embryo_stage")) %>%
  mutate(proportion = if_else(is.na(proportion), 0, proportion))
head(merged_data_updated)

grouped_df <- merged_data_updated %>%
  group_by(embryo_stage, cluster_names) %>%
  summarise(
    proportion = first(proportion), 
    weighted_mean = first(weighted_mean),
    closest_stage = first(closest_stage),
    actual_numeric = first(actual_numeric),
    cluster_germlayer = first(cluster_germlayer)
  ) %>%
  ungroup()  
head(grouped_df)

sum <- grouped_df %>% group_by(embryo_stage) %>% summarize(wm=sum(proportion))

grouped_df <- grouped_df %>%
  arrange(weighted_mean, cluster_germlayer, cluster_names) %>%
  mutate(order_number = (match(embryo_stage, unique(embryo_stage))))
grouped_df$cluster_names <- factor(grouped_df$cluster_names, levels = unique(grouped_df$cluster_names))

# stacked plots
test_stacked_plot <- ggplot(grouped_df, aes(x = factor(weighted_mean), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~closest_stage, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45) )
ggsave("test_stacked_cellstates_ordered_by_prediction_stage_means_LT_10.pdf",width=21, test_stacked_plot)

# confusion matrix.
conf_mat <- confusionMatrix(factor(means$closest_stage), factor(means$actual_numeric))
print(conf_mat)

confMat_melt <- as.data.frame(as.table(conf_mat$table/colSums(conf_mat$table)))
confMat_melt$Prediction <- factor(confMat_melt$Prediction, levels = stages)
confMat_melt$Reference <- factor(confMat_melt$Reference, levels = stages)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Prediction), y =factor(Reference), fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels=named_stages)+
  scale_y_discrete(labels=named_stages)+
  theme_minimal() +
  labs(x = "Prediction", y = "Stage", fill = "Normalized Frequency")
ggsave("LT_10_confusion.pdf",width=10, confMat_plot)

# Classification rate for stages
correct_predictions <- sum(factor(means$actual_numeric) == factor(means$closest_stage))
total_predictions <- nrow(means)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages using weighted means for each embryo (red=",red,", kfilter=",kfilter,"): ", classification_rate))

###################################################################################################################
# means and medians of all prediction scores for each cell and let highest prediction score decide for embryo stage
###################################################################################################################

mean_pred_scores <- data %>% group_by(stage, embryo) %>% 
  summarize(
    WT_90 = mean(prediction.score.WT_90, na.rm = TRUE),
    WT_85 = mean(prediction.score.WT_85, na.rm = TRUE),
    WT_80 = mean(prediction.score.WT_80, na.rm = TRUE),
    WT_75 = mean(prediction.score.WT_75, na.rm = TRUE),
    WT_70 = mean(prediction.score.WT_70, na.rm = TRUE),
    WT_65 = mean(prediction.score.WT_65, na.rm = TRUE)
  )
mean_pred_scores <- mean_pred_scores %>% ungroup() %>%
  mutate(WT_max = pmax(WT_90, WT_85, WT_80, WT_75, WT_70, WT_65, na.rm = TRUE),
  WT_max_col = apply(select(., starts_with("WT_")), 1, function(x) names(which.max(x))) # Find which WT_ column has the max value
  )

# Confusion:
conf_mat <- confusionMatrix(factor(mean_pred_scores$WT_max_col), factor(mean_pred_scores$stage))
print(conf_mat)

# Classification rate for stages
correct_predictions <- sum(factor(mean_pred_scores$WT_max_col) == factor(mean_pred_scores$stage))
total_predictions <- nrow(mean_pred_scores)
classification_rate <- correct_predictions / total_predictions
print(paste("Classification rate for stages using maximal stage prediction means for each embryo (red=",red,", kfilter=",kfilter,"): ", classification_rate))

predicted_numeric <- as.double(paste0(substr(mean_pred_scores$WT_max_col, 4, 4),".",substr(mean_pred_scores$WT_max_col, 5, 5)))
actual_numeric <- as.double(paste0(substr(mean_pred_scores$stage, 4, 4),".",substr(mean_pred_scores$stage, 5, 5) ))

mae <- mae(actual_numeric, predicted_numeric)
mse <- mse(actual_numeric, predicted_numeric)

print(paste0("RMSE = ", sqrt(mse)))

# Confusion Matrix:
confMat_melt <- as.data.frame(as.table(conf_mat$table/colSums(conf_mat$table)))
confMat_melt$Prediction <- factor(confMat_melt$Prediction, levels = stages)
confMat_melt$Reference <- factor(confMat_melt$Reference, levels = stages)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]

# Create the plot
confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Prediction), y =factor(Reference), fill = Freq)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red") +
  scale_x_discrete(labels=named_stages)+
  scale_y_discrete(labels=named_stages)+
  theme_minimal() +
  labs(x = "Prediction", y = "Stage", fill = "Normalized Frequency")
ggsave("LT_10_max_stage_pred_means_confusion.pdf",width=10, confMat_plot)

mean_pred_scores$embryo_stage <- paste0(mean_pred_scores$embryo, "__", mean_pred_scores$stage)
mean_pred_scores <- mean_pred_scores %>% left_join(grouped_df, by="embryo_stage")

sum <- mean_pred_scores %>% group_by(embryo_stage) %>% summarize(wm=sum(proportion))

mean_pred_scores <- mean_pred_scores %>%
  arrange(weighted_mean, cluster_germlayer, cluster_names) %>%
  mutate(order_number = (match(embryo_stage, unique(embryo_stage))))
mean_pred_scores$cluster_names <- factor(mean_pred_scores$cluster_names, levels = unique(mean_pred_scores$cluster_names))
mean_pred_scores$embryo_num <- as.numeric(parse_number( mean_pred_scores$embryo))
mean_pred_scores <- mean_pred_scores %>%
  arrange(cluster_germlayer, cluster_names) 
mean_pred_scores$cluster_names <- factor(mean_pred_scores$cluster_names, levels = unique(mean_pred_scores$cluster_names))
mean_pred_scores$predicted_numeric <- as.double(paste0(substr(mean_pred_scores$WT_max_col, 4, 4),".",substr(mean_pred_scores$WT_max_col, 5, 5)))
mean_pred_scores$actual_numeric  <- as.double(paste0(substr(mean_pred_scores$stage, 4, 4),".",substr(mean_pred_scores$stage, 5, 5)))
mean_pred_scores$embryo_num <- paste0(mean_pred_scores$embryo_num, "_", named_stages[mean_pred_scores$stage])

# stacked plots with WT as reference 
test_stacked_plot <- ggplot(mean_pred_scores, aes(x = factor((embryo_num)), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo Stage Prediction", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~predicted_numeric, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))
ggsave("test_stacked_cellstates_prediction_LT_10.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
mean_pred_scores$delta_pred <- mean_pred_scores$predicted_numeric - mean_pred_scores$actual_numeric
delta_df <- mean_pred_scores %>% group_by(embryo_stage,embryo_num,actual_numeric, predicted_numeric, delta_pred) %>%  summarize()

p <- ggplot(delta_df, aes(x = factor(embryo_num), fill = delta_pred)) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0,limits = c(-2.5, 2.5)) + # Color transition
  labs(x = "", y = "", fill = "Predction-Actual") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lp <- get_legend(p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo Stage Prediction",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())
complete <- plot_grid(test_stacked_plot,p, ncol = 1, rel_heights=c(4,1),align = "h") 
#ggsave(paste0("/project/PRCcomp/images/DevDelay_LT/test_stacked_05delta_complete.pdf"),width=23,complete)

# To align delta and stacked plot with actual (true) heatmap 
a <- ggplot(delta_df, aes(x = factor(embryo_num), fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual Stage", values=color_palette) +
  labs(x = "", y = "", fill = "Actual Stage") +
  theme_minimal() + 
  facet_grid(~predicted_numeric, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 

la <- get_legend(a)
complete <- plot_grid(test_stacked_plot,p,a, ncol = 1, rel_heights=c(4,1,1),align = "h") 
#ggsave(paste0("test_stacked_05delta_actual_complete_LT_10.pdf"),width=23,complete)

delta_df$classification <- ifelse(delta_df$actual_numeric==delta_df$predicted_numeric,"True","False")

b <- ggplot(delta_df, aes(x = factor(embryo_num), fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("False"="red","True"="olivedrab")) +
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
ggsave(paste0("test_stacked_05delta_actual_classification_complete_LT_10.pdf"),width=23,height=9,complete)
ggsave(paste0("test_stacked_05delta_actual_classification_complete_legend.pdf"),ggarrange(lb,lp,la))


#########################################################################################################
# medians of all prediction scores for each cell -> use highest prediction score for stage prediction
#########################################################################################################
mean_pred_scores <- data %>% group_by(stage, embryo) %>% 
  summarize(
    WT_90 = median(prediction.score.WT_90, na.rm = TRUE),
    WT_85 = median(prediction.score.WT_85, na.rm = TRUE),
    WT_80 = median(prediction.score.WT_80, na.rm = TRUE),
    WT_75 = median(prediction.score.WT_75, na.rm = TRUE),
    WT_70 = median(prediction.score.WT_70, na.rm = TRUE),
    WT_65 = median(prediction.score.WT_65, na.rm = TRUE)
  )
mean_pred_scores <- mean_pred_scores %>% ungroup() %>%
  mutate(WT_max = pmax(WT_90, WT_85, WT_80, WT_75, WT_70, WT_65, na.rm = TRUE),
  WT_max_col = apply(select(., starts_with("WT_")), 1, function(x) names(which.max(x))) # Find which WT_ column has the max value
  )
conf_mat <- confusionMatrix(factor(mean_pred_scores$WT_max_col), factor(mean_pred_scores$stage))
print(conf_mat)
# same accuracy with mean()



