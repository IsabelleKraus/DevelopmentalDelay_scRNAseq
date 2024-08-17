library(dplyr)
library(Seurat)
library(MASS)
library(ggplot2)
library(caret)
library(Metrics)
library(viridis)
library(tidyverse)
library(scales) 
library(patchwork)
library(ggpubr)


load("WT.Robj")

# Add annotations to WT object
cluster_order <- c(2, 38, 6, 53, 41, 10, 48, 17, 32, 37, 13, 52, 3, 44, 21, 34, 14, 45, 46, 50, 28, 4, 42, 54, 0, 30, 16, 47, 20, 18, 25, 33, 7, 9, 24, 31, 8, 49, 51, 43, 29, 5, 19, 23, 35, 26, 22, 39, 27, 11, 15, 1, 40, 36, 12)
cluster_names <- c('2_Epiblast', '38_PostProxPrimStreak', '6_EarlyMeso', '53_PGC', '41_BodyWall', '10_ProxExE', '48_DiffExE', '17_DiffExE', '32_Trophoblasts', '37_DefEndo', '13_PrimEndo', '52_YolkSac2', '3_YolkSac1', '44_ParietalEndo', '21_GutEndo', '34_HemogenEndothel', '14_Angioblast', '45_MonocMacrophProg', '46_EMP', '50_HematoEndoProg', '28_PrimBloodProg', '4_Blood_early', '42_PrimFetal_Blood', '54_PrimFetal_Blood', '0_Blood_late', '30_Xmeso_early', '16_Allantois', '47_AmnioticMeso3', '20_AmnioticMeso2', '18_AmnioticMeso1', '25_PresomMeso', '33_EarlyMeso_post', '7_PostLat_IntermMeso', '9_SecHF', '24_PharyArchMeso', '31_PrimHF', '8_Somites', '49_GenitourMeso', '51_Artifact', '43_NodeNotochord', '29_IndEpi_early', '5_IndEpi_late', '19_SurfaceEcto', '23_NeuralRidge_post', '35_NeuralRidge_ant', '26_PrimStreak', '22_NMP', '39_PreneuralPlate_post', '27_PreneuralPlate_ant', '11_Forebrain', '15_Midbrain', '1_NeuralPlate', '40_FloorPlate', '36_MotorNeuron', '12_NeuralCrest')
cluster_colors <- c('#000000', '#474747', '#a0a0a0', '#DB083D', '#722c61', '#fcb7b7', '#f4595c', '#db484a', '#5d0000', '#f6c100', '#f7a800', '#f89e00', '#f98e00', '#fb7100', '#f5f200', '#fad4a6', '#deb88c', '#cba47a', '#bf986f', '#b28b63', '#a57d57', '#98704b', '#8b633f', '#7c5432', '#673f1f', '#2cfeff', '#5bdffe', '#7acafe', '#a9abfd', '#da8afc', '#d767ff', '#a746d0', '#9339bd', '#7927a4', '#62188e', '#3f006c', '#166788', '#099be2', '#492b00', '#1142fc', '#c2ecb3', '#99e2a9', '#4ad094', '#25c88a', '#00bf80', '#00ffda', '#00efca', '#00e1bc', '#00d2ac', '#00c19a', '#00ac83', '#009e75', '#008d62', '#00784b', '#003d0e')
cluster_group <- c('Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Placenta', 'Placenta', 'Placenta', 'Placenta', 'Embryo', 'YolkSac', 'YolkSac', 'YolkSac', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo')
cluster_germlayer <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'Meso', 'XEcto', 'XEcto', 'XEcto', 'XEcto', 'Endo', 'XEndo', 'XEndo', 'XEndo', 'XEndo', 'Endo', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'XEndo', 'Meso', 'Epiblast', 'Epiblast', 'Ecto', 'Ecto', 'Ecto', 'Meso', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto')
germlayer_map <- setNames(cluster_germlayer, cluster_names)

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

# Define a set of distinct colors for each stage (same as in PCA)
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames(stage_colors, e_days)
color_palette <- setNames(stage_colors, c("6.5","7","7.5","8","8.5","9"))
color_mapping2 <- setNames( stage_colors, c("6.5","7.0","7.5","8.0","8.5","9.0"))

# Subset WT data into train and test datasets
testembryos <- c("embryo1", "embryo4", "embryo6")
test <- subset(WT, subset = embryo == testembryos)
train <- subset(WT, subset = embryo != testembryos[1] & embryo != testembryos[2] & embryo != testembryos[3] )  

Idents(train) <- train$stage
all_markers <- list()

for(s in stages) {
  # Find marker genes for cells of a specific stage compared to all others
  markers <- FindMarkers(train, only.pos=TRUE, ident.1 = s, min.pct = 0.05, logfc.threshold = 0.25)
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
length(selected_genes) #should be 905
write.csv(data.frame(Gene = selected_genes), "pos_marker_genes_3embryos.csv", row.names = FALSE)
#selected_genes <- read.csv( "pos_marker_genes_3embryos.csv", header=TRUE)$Gene

# Subset to include only selected genes
train_subset <- subset(train, features = selected_genes)

# Extract data and transpose it
data_for_lda <- GetAssayData(train_subset, slot = "data")
data_for_lda <- as.matrix(data_for_lda)
data_for_lda <- t(data_for_lda)  # transpose

class_labels <- train_subset@meta.data$stage

##############################################################################################################################
# Running LDA (Training)
##############################################################################################################################

lda_result <- lda(data_for_lda, grouping = class_labels, prior=rep(1/6,6))
save(lda_result, file = "lda_result_3embryos.RData")

# Which marker genes are highly LD1-ranked:
sort(abs(lda_result$scaling[,"LD1"]))

lda_prediction <- predict(lda_result, data_for_lda)
str(lda_prediction)

# Extract the LDA scales
lda_scales <- as.data.frame(lda_prediction$x)
# Add the stages to the scores dataframe
lda_scales$stage <- class_labels

proportion_trace <- lda_result$svd^2 / sum(lda_result$svd^2)
print(proportion_trace)
lda_scales$stage <- factor(lda_scales$stage, levels = names(named_stages), labels = named_stages)
stage_labels <- setNames(c("6.5", "7.0", "7.5", "8.0", "8.5", "9.0"), e_days)

# Create the plot for LD1 vs LD2
p1 <- ggplot(lda_scales[sample(1:nrow(lda_scales)), ], aes(x = LD1, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD1 (", round(proportion_trace[1]*100), "%)"), 
       y = paste0("LD2 (", round(proportion_trace[2]*100), "%)"), 
       color = "Embryonic\nstage (E)") + 
  scale_color_manual(values = stage_colors, labels=stage_labels)+
  theme_classic() +
  theme( axis.text = element_text(size = 24, color="black"),  # Increase axis tick text size
    strip.text = element_text(size = 24), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_text(size = 24))+
    guides(color= guide_legend(override.aes = list(size = 3)))
ggsave("LDA_3e_plot.pdf", p1, width = 8, height = 7)

# Create the plot for LD3 vs LD2
p2 <- ggplot(lda_scales[sample(1:nrow(lda_scales)),], aes(x = LD3, y = LD2, color = stage)) + 
  geom_point() +  
  theme_minimal() +  
  theme(legend.position = "right") +  
  labs(x = paste0("LD3 (", round(proportion_trace[3]*100), "%)"), 
       y = paste0("LD2 (", round(proportion_trace[2]*100), "%)"), 
       color = "Stage") + 
  scale_color_manual(values = stage_colors, labels = stage_labels)+
  theme_classic()  +
  theme( axis.text.x = element_text(size = 16, color="black"),  # Increase x axis tick text size
    axis.text.y = element_text(size = 16, color="black"),  # Increase y axis tick text size
    strip.text = element_text(size = 18), axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size = 18))
ggsave("LDA_3e_plot_LD2_LD3.pdf", p2, width = 7, height = 7)

lds<-seq(1,5,1)
# Plotting the proportion of trace as a scree plot
p <- ggplot(data.frame(lds,proportion_trace), aes(x=lds, y=proportion_trace))+
    geom_point(size=3)+
    geom_line()+
    xlab("Number of Linear Discriminants")+
    ylab("Proportion of trace")+
    theme_classic() +
    theme( axis.text.x = element_text(size = 24, color="black"),  # Increase x axis tick text size
    axis.text.y = element_text(size = 24, color="black"),  # Increase y axis tick text size
    strip.text = element_text(size =24), axis.title = element_text(size=24), legend.text = element_text(size=24) ,legend.title = element_text(size = 24))
ggsave("LDA_3e_scree_plot.pdf", p, width=7,height=7)

#########################################################################################################################################################
# TEST (Validation)
#########################################################################################################################################################
test_subset <- subset(test, features = selected_genes)

# Extract data
test_for_lda <- t(as.matrix(GetAssayData(test_subset, slot = "data")))
class_labels <- test_subset@meta.data$stage
test_lda_prediction <- predict(lda_result, test_for_lda)
str(test_lda_prediction)
# Extract the scales
test_lda_scales <- as.data.frame(test_lda_prediction$x)
# Add the stages to the scores dataframe
test_lda_scales$stage <- class_labels
test_lda_scales$stage <- factor(test_lda_scales$stage, levels = names(named_stages), labels = named_stages)
test_subset$pred <- test_lda_prediction$class
conf_mat <- confusionMatrix(factor(test_subset$pred), factor(test_subset$stage))
print(conf_mat)

# Confusion Matrix:
confMat_melt <- as.data.frame(sweep(as.table(conf_mat$table), 2, colSums(conf_mat$table), FUN="/"))
confMat_melt$Prediction <- factor(confMat_melt$Prediction)
confMat_melt$Reference <- factor(confMat_melt$Reference)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]
named_stages2 <- c("WT_65" = "6.5","WT_70" = "7.0", "WT_75" = "7.5","WT_80" = "8.0", "WT_85" = "8.5","WT_90" = "9.0")

# confusion matrix plot
confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Reference), y =factor(Prediction, levels = rev(levels(factor(Prediction)))), fill = Freq)) +
  geom_tile(color="grey") +
  scale_fill_gradient(low = "white", high = "navy") +
  scale_x_discrete(labels=named_stages2)+
  scale_y_discrete(labels=rev(named_stages2))+
  theme_minimal() +
  labs(x = "Embryonic stage (E)", y = "Cell prediction (E)", fill = "Normalized \nfrequency")+
  theme(axis.text=element_text(size = 24, color="black"),  
        strip.text = element_text(size = 24), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_text(size = 24, angle=90)) 
ggsave("LDA_3e_Confmat_cells_normalized.pdf",width=9, confMat_plot)

test_subset$predicted_numeric <- as.double(paste0(substr(test_subset$pred, 4, 4),".",substr(test_subset$pred, 5, 5)))
test_subset$actual_numeric  <- as.double(paste0(substr(test_subset$stage, 4, 4),".",substr(test_subset$stage, 5, 5)))
test_subset$delta_pred <- test_subset$predicted_numeric - test_subset$actual_numeric

# Root mean square error
mse <- mse(test_subset$actual_numeric, test_subset$predicted_numeric)
print(paste0("RMSE = ", sqrt(mse)))

test_subset$e_day <- named_stages[test_subset$stage]
test_subset$pred_e <- named_stages[test_subset$pred]
test_subset$pred_d <- named_stages2[test_subset$pred]
train$e_day <- named_stages[train$stage]
train$day <- named_stages2[train$stage]

# UMAP for Train Embryos 
umap_train <- DimPlot(train, cols=color_mapping2, group.by="day", pt.size=1, shuffle=TRUE, seed=42)+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(), 
    legend.text = element_text(size=24),
    legend.title = element_text(size=24)
  ) +
  guides(color = guide_legend(title = "Embryonic\nstage (E)",override.aes = list(size = 3))) +
  ggtitle("")
ggsave("LDA_3e_umap_train.pdf",width=8, umap_train)

# UMAP for Test Embryos
umap_pred <- DimPlot(test_subset, cols=color_mapping2, pt.size=1, shuffle=TRUE, seed=42, group.by="pred_d")+
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size=24),
    legend.title = element_text(size=24)) +
    guides(color = guide_legend(title = "Prediction (E)",override.aes = list(size = 3))) +
    ggtitle("")
ggsave("LDA_3e_umap_test_pred.pdf",width=8, umap_pred)

# UMAP for delta values 
delta_colors <- c(
  "0" = "gray",        
  "-1.5" = "royalblue",    
  "-1" = "cornflowerblue", 
  "-0.5" = "lightskyblue",              
  "0.5" = "lightcoral",    
  "1" = "indianred",       
  "1.5" = "darkred"        
)
alpha_value <- 0.3  # Set the alpha value (0 = fully transparent, 1 = fully opaque)
delta_colors[1] <- alpha(delta_colors[1], alpha_value)

umap_delta <- DimPlot(test_subset, cols = delta_colors, pt.size=2,  shuffle=TRUE, seed=42,group.by="delta_pred")+theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size=24),
    legend.title = element_text(size=24, angle=90)) +
    guides(color = guide_legend(title = "Prediction error\n(Prediction - Actual stage)",override.aes = list(size = 3)))+
  ggtitle("")
ggsave("LDA_3e_umap_delta_predictions.pdf",width=8,umap_delta)

# Which cell states were misclassified:
delta_test <- test_subset@meta.data %>% filter(delta_pred!=0)
table(delta_test$delta_pred, delta_test$cluster_name)
sort(colSums(table(delta_test$delta_pred, delta_test$cluster_name)))

# Which embryos were misclassified:
sort(colSums(table(delta_test$delta_pred, paste(delta_test$embryo,delta_test$stage))))
embryo_nums <- test_subset@meta.data %>%
  group_by(paste(embryo,stage)) %>%
  summarise(n = n()) 
  
sort(colSums(table(delta_test$delta_pred, paste(delta_test$embryo,delta_test$stage)))/embryo_nums$n)

# Posterior probabilities give a measure of certainty about each prediction:
# Find amount of predictions with high <0.9 and low <=0.9 certainty
high_certainty_filter <- apply(test_lda_prediction$posterior, 1, max) > 0.9
print(table(high_certainty_filter))
test_subset <- AddMetaData(test_subset, metadata = as.data.frame(high_certainty_filter))
test_subset@meta.data$classification <- ifelse(test_subset@meta.data$pred  == test_subset@meta.data$stage, "TRUE", "FALSE")

# Calculate the intersections
certain <- sum(test_subset@meta.data$high_certainty_filter == TRUE)
uncertain <- sum(test_subset@meta.data$high_certainty_filter == FALSE)
both_true <- sum(test_subset@meta.data$classification == TRUE & test_subset@meta.data$high_certainty_filter == TRUE)
class_true_filter_false <- sum(test_subset@meta.data$classification == TRUE & test_subset@meta.data$high_certainty_filter == FALSE)
class_false_filter_true <- sum(test_subset@meta.data$classification == FALSE & test_subset@meta.data$high_certainty_filter == TRUE)
both_false <- sum(test_subset@meta.data$classification == FALSE & test_subset@meta.data$high_certainty_filter == FALSE)

input <- data.frame(x = c("True","True", "False", "False"), 
  vals=c(both_true, class_true_filter_false, class_false_filter_true, both_false), 
  filter=c("p > 0.9", "p <= 0.9","p > 0.9", "p <= 0.9"))
input <- input %>% mutate(percentage= vals / sum(vals))

# stacked bar plot + analysis of classification  
st <- ggplot(input, aes(x = x, y = percentage*100, fill = filter)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Classification",
       y = "Percentage of cells",
       fill = "Posterior probability") +
  theme_minimal()+
  theme(axis.text.x=element_text(size = 12),    
        axis.text.y = element_text(size = 12),  
        strip.text = element_text(size = 14), axis.title = element_text(size=16) )
ggsave("LDA_3e_Classification_stacked.pdf", width=7, plot = st)

# Distribution of maximal Posterior values 
max_posterior <- apply(test_lda_prediction$posterior, 1, max)
test_subset <- AddMetaData(test_subset, metadata = as.data.frame(max_posterior))

box <- ggplot(test_subset@meta.data, aes(x=factor(classification, levels = rev(levels(factor(classification)))), y=max_posterior))+
    geom_boxplot(fill="lightgrey",  outlier.shape=1)+
    labs(x= "Classification", y="Posterior probability")+
    scale_x_discrete(name = "Classification", labels = c("TRUE"="True", "FALSE"="False")) +
    theme_minimal()+
    theme_classic() +
    theme( axis.text.x = element_text(size = 24, color="black"),  # Increase x axis tick text size
    axis.text.y = element_text(size = 24, color="black"),  # Increase y axis tick text size
    strip.text = element_text(size = 24), axis.title = element_text(size=24), legend.text = element_text(size=24), legend.title = element_text(size = 24))
ggsave("LDA_3e_Classification_boxplot.pdf", width=4,box)

# Classification of validation embryos 
test_subset <- AddMetaData(test_subset, metadata = as.data.frame(test_lda_prediction$posterior))
stages_numeric <- c(6.5, 7.0, 7.5, 8.0, 8.5, 9.0)

# weighted sum of stages
weighted_sums <- apply(test_lda_prediction$posterior, 1, function(x) sum(x * stages_numeric))
test_subset <- AddMetaData(test_subset, metadata = as.data.frame(weighted_sums))

# predicted_numeric stage for each embryo:
means <- test_subset@meta.data %>% group_by(actual_numeric, embryo) %>%  summarise(mean = mean(weighted_sums)) %>%
    ungroup()
means$closest_stage <- sapply(means$mean, function(x) {
  stages_numeric[which.min(abs(x - stages_numeric))]
})

# confusion matrix + statistics
conf_mat <- confusionMatrix(factor(means$closest_stage), factor(means$actual_numeric))
print(conf_mat)

# Confusion Matrix:
confMat_melt <- as.data.frame(as.table(conf_mat$table))
confMat_melt$Prediction <- factor(confMat_melt$Prediction)
confMat_melt$Reference <- factor(confMat_melt$Reference)
confMat_melt <- confMat_melt[order(confMat_melt$Reference, confMat_melt$Prediction),]
named_stages2 <- c("6.5" = "E6.5","7" = "E7.0", "7.5" = "E7.5","8" = "E8.0", "8.5" = "E8.5","9" = "E9.0")

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Reference), y =factor(Prediction, levels=rev(levels(factor(Prediction)))), fill = Freq)) +
  geom_tile(color="grey") +
  scale_fill_gradient(low = "white", high = "navy") +
  scale_x_discrete(labels=named_stages2)+
  scale_y_discrete(labels=rev(named_stages2))+
  theme_minimal() +
  labs(x = "Embryonic stage (E)", y = "Embryo prediction (E)", fill = "Frequency")+
  theme(
    axis.text = element_text(size = 24, color="black"),  
    strip.text = element_text(size = 24), axis.title = element_text(size=24),
    legend.text = element_text(size=24), legend.title = element_text(size = 24, angle=90))
ggsave("LDA_3e_Confmat_embryos.pdf",width=9, confMat_plot)

confMat_plot <- ggplot(data = confMat_melt, aes(x = factor(Reference), y =factor(Prediction, levels=rev(levels(factor(Prediction)))), fill = Freq/3)) +
  geom_tile(color="grey") +
  scale_fill_gradient(low = "white", high = "navy") +
  scale_x_discrete(labels=named_stages2)+
  scale_y_discrete(labels=rev(named_stages2))+
  theme_minimal() +
  labs(x = "Embryonic stage (E)", y = "Embryo prediction (E)", fill = "Normalized\nFrequency")+
  theme(
    axis.text = element_text(size = 24, color="black"),  
    strip.text = element_text(size = 24), axis.title = element_text(size=24),
    legend.text = element_text(size=24), legend.title = element_text(size = 24, angle=90))
ggsave("LDA_3e_Confmat_embryos_normalized.pdf",width=9, confMat_plot)

# RMSE of continuous stages
mse <- mse(means$actual_numeric, means$mean)
print(paste0("RMSE = ", sqrt(mse)))

# RMSE of discrete stages
mse <- mse(means$actual_numeric, means$closest_stage)
print(paste0("RMSE = ", sqrt(mse)))

test_subset@meta.data$embryo_name <- paste(test_subset@meta.data$embryo, test_subset@meta.data$actual_numeric, sep="_")

# cell state poportions 
all_combinations <- expand.grid(embryo = test_subset@meta.data$embryo_name, cluster_names = cluster_names)
total_cells <- test_subset@meta.data %>% group_by(embryo_name) %>%  summarise(total = n(), .groups = 'drop')
proportions <- test_subset@meta.data %>%
  left_join(total_cells, by = "embryo_name") %>%
  group_by(cluster_names, embryo_name, total, sex) %>% # Include 'total' in group_by
  summarise(count = n(), .groups = 'drop') %>%
  mutate(proportion = count / total)
proportions_complete <- all_combinations %>%
  left_join(proportions, by = c("cluster_names", "embryo"="embryo_name")) %>%
  replace_na(list(count = 0, proportion = 0)) 

# Stacked Plot for each test embryo: 
proportions_complete <- proportions_complete %>%
  mutate(cluster_germlayer = germlayer_map[cluster_names]) 
  arrange(cluster_germlayer, cluster_names) %>%
  mutate(embryo_num = parse_number(embryo))
proportions_complete$cluster_names <- factor(proportions_complete$cluster_names, levels = unique(proportions_complete$cluster_names))

means$embryo_name <- paste(means$embryo, means$actual_numeric, sep="_")
merged_data <- proportions %>%
  left_join(means %>% dplyr::select(embryo_name,closest_stage,mean), 
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
  group_by(embryo_name, mean, closest_stage) %>%  # Group by both order_number and stage
  summarize(sex_label = substr(first(sex),1,1), .groups = 'drop')  # Get the first sex label for each group
embryo_sex_label$embryo_name <- factor(embryo_sex_label$embryo_name, levels = unique(embryo_sex_label$embryo_name))

test_stacked_plot <- ggplot(merged_data, aes(x = factor(round(mean,4)), y = proportion*100, fill = cluster_names)) +
        geom_bar(stat = "identity", position = "stack") +
        #facet_grid(scale="free", space="free_x") +
        theme_minimal() +
        scale_x_discrete()+
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryonic stage prediction", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        facet_grid(~closest_stage, scales="free",space="free_x")+
        #ggtitle("E8.5 WT")+
        theme(legend.position = "none", axis.text.x=element_text(angle=45, hjust=1))+
        theme(axis.text=element_text(size = 30, color="black"),  # Increase y axis tick text size
        strip.text = element_text(size = 30), axis.title = element_text(size=30), legend.text = element_text(size=30), legend.title = element_text(size = 30)) 
ggsave("LDA_3e_test_stacked_cellstates_ordered_by_pred_stage_means_uniform.pdf",width=21, test_stacked_plot)

test_stacked_plot <- test_stacked_plot + labs(x="Embryonic stage prediction",y="") + theme( axis.text.y=element_blank(), 
      axis.ticks.y=element_blank(),axis.title.x=element_text(size=30), strip.text.x = element_blank())

delta_df <- merged_data %>% group_by(embryo_name,closest_stage, mean) %>%  summarize()
delta_df$classification <- "True"

b <- ggplot(delta_df, aes(x = factor(round(mean,4)), fill = factor(classification))) +
  geom_bar(stat = "identity", aes(y = 1)) + 
  scale_fill_manual(name="Classification", values=c("False"="red","True"="olivedrab")) +
  labs(x = "", y = "", fill = "Classification") +
  theme_minimal() + 
  facet_grid(~closest_stage, scales="free",space="free_x") +
  guides(color = guide_legend(title = "Classification",override.aes = list(size = 5))) +
  theme(axis.text.y=element_blank(),  legend.text = element_text(size=24),legend.title = element_text(size = 24),  strip.text = element_text(size = 30),
      axis.ticks.y=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm")) 
b <- b + rremove("xlab") + rremove("legend")

complete <- b/test_stacked_plot + plot_layout(heights = c(0.2,4))
ggsave(paste0("LDA_3e_test_stacked_delta_actual_classification_complete.pdf"),width=23,height=9,complete)

