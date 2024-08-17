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
library(Matrix)  
library(MASS)
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

# Define a set of distinct colors for each stage 
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames(stage_colors, e_days)
color_palette <- setNames(stage_colors, c("6.5","7","7.5","8","8.5","9"))

red_palette <- c(
  "6.5" = "#FFCCCC",  # Lighter red
  "7" = "#FF9999",
  "7.5" = "#FF6666",
  "8" = "#FF3333",
  "8.5" = "#FF0000",
  "9" = "#CC0000"   # Darker red
)

# Determine subsets (10% validation, 90% train)
set.seed(42) # Setting seed for reproducibility
cells_90 <- sample(all_cells, size = length(all_cells) * 0.9)

# Split the WT object into two
WT_90 <- subset(WT, cells = cells_90)
WT_10 <- subset(WT, cells = setdiff(all_cells, cells_90))

###################################################################################################
# Training
###################################################################################################

Idents(WT_90) <- WT_90$stage
all_markers <- list()

for(s in stages) {
  # Find marker genes for cells of a specific stage compared to all others
  markers <- FindMarkers(WT_90, ident.1 = s,only.pos=TRUE, min.pct = 0.05, logfc.threshold = 0.25)
  all_markers[[s]] <- markers
}

# Filter markers based on p-value threshold (p < 0.05) for each stage
all_markers <- lapply(all_markers, function(markers) {
  markers[markers$p_val < 0.05, ]
})

# Check entries in the list for each stage
lapply(all_markers, head)

# Create a vector of unique markers from all stages
selected_genes <- unique(unlist(lapply(all_markers, rownames)))
write.csv(data.frame(Gene = selected_genes), "pos_marker_genes_10.csv", row.names = FALSE)

# Subset to include only selected genes
WT_90_subset <- subset(WT_90, features = selected_genes)

# Extract data 
data_for_lda <- GetAssayData(WT_90_subset, slot = "data")
data_for_lda <- t(data_for_lda) 
data_for_lda <- as.matrix(data_for_lda)
class_labels <- WT_90_subset@meta.data$stage

# Running LDA
lda_result <- lda(data_for_lda, grouping = class_labels, prior=rep(1/6,6))
save(lda_result, file = "lda_result_10.RData")

lda_prediction <- predict(lda_result, data_for_lda)
str(lda_prediction)

# Extract the scores
lda_scores <- as.data.frame(lda_prediction$x)

# Add the stages to the scores dataframe
lda_scores$stage <- class_labels

# proportion of trace: reflects proportion of the total discriminative information captured by each LD
proportion_trace <- lda_result$svd^2 / sum(lda_result$svd^2)
print(proportion_trace)

lda_scores$stage <- factor(lda_scores$stage, levels = names(named_stages), labels = named_stages)

# plot for LD1 vs LD2
p1 <- ggplot(lda_scores[sample(1:nrow(lda_scores)), ], aes(x = LD1, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD1 (", round(proportion_trace[1]*100), "%)"), 
       y = paste0("LD2 (", round(proportion_trace[2]*100), "%)"), 
       color = "Stage") + 
  scale_color_manual(values = stage_colors)
ggsave("LDA_10_plot.pdf", p1, width = 7, height = 7)

# plot for LD3 vs LD2
p2 <- ggplot(lda_scores[sample(1:nrow(lda_scores)),], aes(x = LD3, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD3 (", round(proportion_trace[3]*100), "%)"), 
       y = paste0("LD2 (", round(proportion_trace[2]*100), "%)"), 
       color = "Stage") + 
  scale_color_manual(values = stage_colors)
ggsave("LDA_10_plot_LD2_LD3.pdf", p2, width = 7, height = 7)

lds<-seq(1,5,1)
# Plotting the proportion of trace as a scree plot
p <- ggplot(data.frame(lds,proportion_trace), aes(x=lds, y=proportion_trace))+
    geom_point()+
    geom_line()+
    xlab("Number of Linear Discriminants")+
    ylab("Proportion of Trace")+
    theme_minimal()
ggsave("LDA_10_scree_plot.pdf", p)

###############################################################################################################
# Validation
###############################################################################################################

WT_10_subset <- subset(WT_10, features = selected_genes)

# Extract data 
test_for_lda <- t(GetAssayData(WT_10_subset, slot = "data"))
test_for_lda <- as.matrix(test_for_lda)
class_labels <- WT_10_subset@meta.data$stage

test_lda_prediction <- predict(lda_result, test_for_lda)
str(test_lda_prediction)

# Extract the scores
test_lda_scores <- as.data.frame(test_lda_prediction$x)
test_lda_scores$stage <- class_labels
test_lda_scores$stage <- factor(test_lda_scores$stage, levels = names(named_stages), labels = named_stages)

# Validation LDA plots
p <- ggplot(test_lda_scores[sample(1:nrow(test_lda_scores)),], aes(x = LD1, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD 1 (", round(proportion_trace[1]*100,1),"%)"), y = paste0("LD2 (", round(proportion_trace[2]*100,1),"%)"), color = "Stage") + 
  scale_color_manual(values = stage_colors)  
ggsave("10_test_LDA_plot.pdf", p, width = 7, height = 7)

p <- ggplot(test_lda_scores[sample(1:nrow(test_lda_scores)),], aes(x = LD3, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD 3 (", round(proportion_trace[3]*100,1),"%)"), y = paste0("LD2 (", round(proportion_trace[2]*100,1),"%)"), color = "Stage") + 
   scale_color_manual(values = stage_colors)  
ggsave("10_test_LDA_plot_LD2_LD3.pdf", p, width = 7, height = 7)

WT_10_subset$pred <- test_lda_prediction$class

# Classification statistics
conf_mat <- confusionMatrix(factor(WT_10_subset$pred), factor(WT_10_subset$stage))
print(conf_mat)

WT_10_subset$predicted_numeric <- as.double(paste0(substr(WT_10_subset$pred, 4, 4),".",substr(WT_10_subset$pred, 5, 5)))
WT_10_subset$actual_numeric  <- as.double(paste0(substr(WT_10_subset$stage, 4, 4),".",substr(WT_10_subset$stage, 5, 5)))
WT_10_subset$delta_pred <- WT_10_subset$predicted_numeric - WT_10_subset$actual_numeric

mae <- mae(WT_10_subset$actual_numeric, WT_10_subset$predicted_numeric)
mse <- mse(WT_10_subset$actual_numeric, WT_10_subset$predicted_numeric)
print(paste0("RMSE = ", sqrt(mse)))

WT_10_subset$e_day <- named_stages[WT_10_subset$stage]
WT_10_subset$pred_e <- named_stages[WT_10_subset$pred]

# UMAPs
# Predicted and actual stages of validation set
umap_pred <- DimPlot(WT_10_subset, cols=color_mapping, group.by="pred_e")
ggsave("LDA_10_umap_lda_pred.pdf",width=8, umap_pred)
umap_real <- DimPlot(WT_10_subset, cols=color_mapping, group.by="e_day")
ggsave("LDA_10_umap_actual_stages.pdf",width=8, umap_real)

# correct and false predicted cells
correct <- WT_10_subset[,WT_10_subset@meta.data$stage == WT_10_subset@meta.data$pred]
false <- WT_10_subset[,WT_10_subset@meta.data$stage != WT_10_subset@meta.data$pred]
umap_correct <- DimPlot(correct, cols=color_mapping, group.by="e_day")
umap_false <- DimPlot(false, cols=color_mapping, group.by="pred_e")
ggsave("LDA_10_umap_correct_predictions.pdf",width=8,umap_correct)
ggsave("LDA_10_umap_false_predictions.pdf",width=8, umap_false)

myColorPalette <- colorRampPalette(c("blue","gray89", "red"))
colors <- myColorPalette(9)
# colors for delta values
delta_colors <- setNames(colors, seq(-2.0,2.0,0.5))

# cells showing prediction error (deltas)
umap_delta <- DimPlot(false,cols=delta_colors,  group.by="delta_pred")
ggsave("LDA_10_umap_delta_false_predictions.pdf",width=8,umap_delta)
umap_delta_all <- DimPlot(WT_10_subset,cols=delta_colors, order=c(-1, 1,-0.5,0.5,0),  group.by="delta_pred")
ggsave("LDA_10_umap_delta_predictions.pdf",width=8,umap_delta_all)


# Posterior probabilities 

# Find amount of predictions with high <0.9 and low <=0.9 certainty
high_certainty_filter <- apply(test_lda_prediction$posterior, 1, max) > 0.9
print(table(high_certainty_filter))
WT_10_subset <- AddMetaData(WT_10_subset, metadata = as.data.frame(high_certainty_filter))
WT_10_subset@meta.data$classification <- ifelse(WT_10_subset@meta.data$pred  == WT_10_subset@meta.data$stage, "TRUE", "FALSE")

# Calculate intersections
certain <- sum(WT_10_subset@meta.data$high_certainty_filter == TRUE)
uncertain <- sum(WT_10_subset@meta.data$high_certainty_filter == FALSE)
both_true <- sum(WT_10_subset@meta.data$classification == TRUE & WT_10_subset@meta.data$high_certainty_filter == TRUE)
class_true_filter_false <- sum(WT_10_subset@meta.data$classification == TRUE & WT_10_subset@meta.data$high_certainty_filter == FALSE)
class_false_filter_true <- sum(WT_10_subset@meta.data$classification == FALSE & WT_10_subset@meta.data$high_certainty_filter == TRUE)
both_false <- sum(WT_10_subset@meta.data$classification == FALSE & WT_10_subset@meta.data$high_certainty_filter == FALSE)

input <- data.frame(x = c("True","True", "False", "False"), 
  vals=c(both_true, class_true_filter_false, class_false_filter_true, both_false), 
  filter=c("p > 0.9", "p <= 0.9","p > 0.9", "p <= 0.9"))
input <- input %>% mutate(percentage= vals / sum(vals))

# stacked bar plot of classification result and posterior values 
st <- ggplot(input, aes(x = x, y = percentage*100, fill = filter)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Classification",
       y = "Percentage of cells",
       fill = "Posterior probability") +
  theme_minimal()
ggsave("LDA_10_Classification_stacked.pdf", width=5, plot = st)

# Distribution of maximal Posterior values 
max_posterior <- apply(test_lda_prediction$posterior, 1, max)
WT_10_subset <- AddMetaData(WT_10_subset, metadata = as.data.frame(max_posterior))

box <- ggplot(WT_10_subset@meta.data, aes(x=classification, y=max_posterior))+
    geom_boxplot()+
    labs(x= "Classification", y="Posterior probability")+
    scale_x_discrete(name = "Classification", labels = c("TRUE"="True", "FALSE"="False")) +
    theme_minimal()
ggsave("LDA_10_Classification_boxplot.pdf",box)

WT_10_subset <- AddMetaData(WT_10_subset, metadata = as.data.frame(test_lda_prediction$posterior))
stages_numeric <- rev(c(9.0,8.5,8.0,7.5,7.0,6.5))

# Compute weighted sum of stages
weighted_sums <- apply(test_lda_prediction$posterior, 1, function(x) sum(x * stages_numeric))
WT_10_subset <- AddMetaData(WT_10_subset, metadata = as.data.frame(weighted_sums))

# predicted_numeric stage for every embryo:
means <- WT_10_subset@meta.data %>% group_by(actual_numeric, embryo) %>%  summarise(mean = mean(weighted_sums)) %>%
    ungroup()
means$closest_stage <- sapply(means$mean, function(x) {
  stages_numeric[which.min(abs(x - stages_numeric))]
})

# Embryo classification statistics
conf_mat <- confusionMatrix(factor(means$closest_stage), factor(means$actual_numeric))
print(conf_mat)

# Confusion Matrix:
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
ggsave("LDA_10_all_pred_scores_confusion.pdf",width=10, confMat_plot)

mae <- mae(means$actual_numeric, means$mean)
mse <- mse(means$actual_numeric, means$mean)
print(paste0("RMSE (using predicted mean values) = ", sqrt(mse)))

mae <- mae(means$actual_numeric, means$closest_stage)
mse <- mse(means$actual_numeric, means$closest_stage)
print(paste0("RMSE (using predicted continuous stages) = ", sqrt(mse)))

WT_10_subset@meta.data <- WT_10_subset@meta.data %>%
  mutate(embryo_stage = paste(embryo, stage, sep = "__"))
WT_10_subset@meta.data <- WT_10_subset@meta.data %>% left_join(means, by=c("actual_numeric", "embryo"))

test_all_combinations <- expand.grid(embryo_stage = unique(WT_10_subset@meta.data$embryo_stage), 
                                cluster_names = cluster_names)
test_total_cells <- WT_10_subset@meta.data %>% group_by(embryo_stage) %>%  summarise(total = n(), .groups = 'drop')

# Calculate proportions for each combination of cluster_name, stage, and embryo
test_proportions <-  WT_10_subset@meta.data %>%
  left_join(test_total_cells, by = c( "embryo_stage")) %>%
  group_by(cluster_names, stage, embryo, embryo_stage, total ) %>% 
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)

merged_data_updated <- WT_10_subset@meta.data %>%
  left_join(test_proportions %>% dplyr::select(cluster_names, embryo_stage, proportion), 
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
    cluster_germlayer = first(cluster_germlayer),
    delta = first(delta_pred)
  ) %>%
  ungroup()  
head(grouped_df)

# just a control (should be 1): 
#sum <- grouped_df %>% group_by(embryo_stage) %>% summarize(wm=sum(proportion))

grouped_df <- grouped_df %>%
  arrange(cluster_germlayer, cluster_names) 
grouped_df$cluster_names <- factor(grouped_df$cluster_names, levels = unique(grouped_df$cluster_names))

# stacked plots of validation embryos
test_stacked_plot <- ggplot(grouped_df, aes(x = factor(round(mean,7)), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo Stage Prediction", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~closest_stage, scales="free",space="free_x")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))
ggsave("LDA_10_test_stacked_cellstates_ordered_by_prediction_stage_means.pdf",width=21, test_stacked_plot)

# To align stacked plot with delta plot: 
grouped_df$delta <- grouped_df$mean - grouped_df$actual_numeric
delta_df <- grouped_df %>% group_by(embryo_stage,actual_numeric,closest_stage, mean, delta) %>%  summarize()

p <- ggplot(delta_df, aes(x = factor(round(mean,7)), fill = delta)) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_gradient2(low = "blue", mid = "gray89", high = "red", midpoint = 0,  limits = c(-2.5, 2.5)) + # Color transition
  labs(x = "", y = "", fill = "Predction-Actual") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),strip.text.x = element_blank(),axis.text.x=element_blank(),plot.margin = unit(c(0, 0, 0, 0), "cm")) 
lp <- get_legend(p)
ggsave("LDA_10_prediction_delta.pdf",width=21,height=2, p)

test_stacked_plot <- test_stacked_plot + labs(x="Embryo Stage Prediction",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(), strip.text.x = element_blank())

# To align delta and stacked plot with actual stages 
a <- ggplot(delta_df, aes(x = factor(round(mean,7)), fill = factor(actual_numeric))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Actual Stage", values=color_palette) +
  labs(x = "", y = "", fill = "Actual Stage") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  theme(axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(), strip.text.x = element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
la <- get_legend(a)

# To align stacked plot with classification result
delta_df$classification <- ifelse(delta_df$actual_numeric==delta_df$closest_stage,"True","False")

b <- ggplot(delta_df, aes(x = factor(round(mean,7)), fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("False"="red","True"="olivedrab")) +
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
ggsave(paste0("LDA_10_test_stacked_delta_actual_classification_complete.pdf"),width=23,height=9,complete)
ggsave(paste0("LDA_10_test_stacked_delta_actual_classification_complete_legend.pdf"),ggarrange(lb,lp,la))

#######################################################################################################################################################################
# 2nd classification: mean over all posteriors per stage (for each embryo) + argmax
#######################################################################################################################################################################

mean_pred_scores <- WT_10_subset@meta.data %>% group_by(stage, embryo) %>% 
  summarize(
    WT_90 = mean(WT_90, na.rm = TRUE),
    WT_85 = mean(WT_85, na.rm = TRUE),
    WT_80 = mean(WT_80, na.rm = TRUE),
    WT_75 = mean(WT_75, na.rm = TRUE),
    WT_70 = mean(WT_70, na.rm = TRUE),
    WT_65 = mean(WT_65, na.rm = TRUE)
  )

# function to find the column name of the maximum value
find_max_col <- function(row) {
  vals <- c(WT_90 = row$WT_90, WT_85 = row$WT_85, WT_80 = row$WT_80, WT_75 = row$WT_75, WT_70 = row$WT_70, WT_65 = row$WT_65)
  # Return name of the maximum value
  names(which.max(vals))
}

mean_pred_scores <- mean_pred_scores %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    WT_max = max(WT_90, WT_85, WT_80, WT_75, WT_70, WT_65, na.rm = TRUE),  # Find maximum value
    WT_max_col = find_max_col(cur_data())  # Use custom function to find the column name of the max value
  ) %>%
  ungroup() 

# Classification statistics
conf_mat <- confusionMatrix(factor(mean_pred_scores$WT_max_col), factor(mean_pred_scores$stage))
print(conf_mat)







