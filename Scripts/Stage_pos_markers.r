library(dplyr)
library(Seurat)
library(tidyverse)
library(ComplexUpset)
library(ggplot2)

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
all_cells <- colnames(WT@assays$RNA@counts)
unique_stages <- unique(WT@meta.data$stage)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, unique_stages)

train <- WT
Idents(train) <- train$stage
all_markers <- list()

for(s in stages) {
  # Find marker genes for cells of a specific stage compared to all others (only positive ones)
  markers <- FindMarkers(train, only.pos=TRUE, ident.1 = s, min.pct = 0.05, logfc.threshold = 0.25)
  all_markers[[s]] <- markers
}

# Check entries in the list for each stage
lapply(all_markers, head)

# Filter markers based on p-value threshold (p < 0.05) for each stage
all_markers <- lapply(all_markers, function(markers) {
  markers[markers$p_val < 0.05, ]
})

# Export individual stage gene lists to CSV
lapply(names(all_markers), function(set_name) {
  write.csv(all_markers[[set_name]], paste0(set_name, "_pos_markers.csv"), row.names = TRUE)
})

# Combine all genes into one data frame with their set names (stage)
gene_lists <- lapply(names(all_markers), function(set_name) {
  data.frame(gene = rownames(all_markers[[set_name]]), set = set_name)
})

combined_genes <- do.call(rbind, gene_lists)

# Count occurrences of each gene across sets
gene_counts <- combined_genes %>% 
  group_by(gene) %>%
  summarize(count = n(), .groups = 'drop')

# Sort the genes by their count of appearances in descending order
sorted_genes <- arrange(gene_counts, desc(count))
print(sorted_genes)

all_markers <- list()

# Import individual stage gene lists from CSV, only the first column
all_markers <- setNames(lapply(stages, function(set_name) {
  # first column contains name of genes 
  df <- read.csv(paste0(set_name, "_pos_markers.csv"),
                 header = TRUE, colClasses = c("character", rep("NULL", 5)))
   df[, 1] <- toupper(df[, 1])
  colnames(df) <- c("Gene")
  df
}), stages)


# 'all_markers' is a list containing various entries including those prefixed with 'WT_'
# export each gene list to a CSV file
lapply(stages, function(set_name) {
  # Extract gene names
  gene_data <- all_markers[[set_name]]$"Gene"   
  # Create a gene name txt file for each set
  filename <- paste0(set_name, "_pos_genes.txt")  # Specify your path
  writeLines(gene_data, filename)
  cat(set_name, "genes exported to", filename, "\n")
})

####################################################################################################
# Generate file names + read
file_names <- sprintf("WT_%d_pos_genes.txt", seq(65, 90, by = 5))

gene_lists <- lapply(file_names, function(file) {
  read.table(file, header = FALSE, stringsAsFactors = FALSE, col.names = c("gene"))
})

# adding column indicating the stage each gene came from
gene_lists <- lapply(seq_along(gene_lists), function(i) {
  mutate(gene_lists[[i]], stage = stages[i])
})
all_genes <- bind_rows(gene_lists)
head(all_genes)

# format suitable for UpSet plot
gene_presence <- all_genes %>%
  mutate(presence = TRUE) %>%
  pivot_wider(names_from = stage, values_from = presence, values_fill = list(presence = FALSE))
print(gene_presence)

gene_presence_b <- gene_presence 
colnames(gene_presence_b)[2:7] <- e_days

# Generate UpSet plot showing stage intersections of marker genes  
up <- upset(gene_presence_b, rev(e_days), name="Stage", base_annotations=list('Intersection size'=intersection_size(counts=FALSE)),sort_sets = FALSE, themes=upset_default_themes(text=element_text(size=24)))
ggsave("Upset_upregulated.pdf", up, width=14)

# Generate gene lists ranked by LD1 scale (absolute value) 
load("lda_result_WT.RData")
lda_scaling <- lda_result$scaling

# CSV files of unique-stage markers for GO enrichment analysis (ranked by LD1 values)
  for (stage in stages) {
    # Filter marker genes that only occur in one stage
    stage_genes <- gene_presence$gene[gene_presence[[stage]] & rowSums(gene_presence[stages]) == 1]
    ld1_values <- lda_scaling[match(stage_genes, toupper(rownames(lda_scaling))),"LD1"]
    print(ld1_values)

    result_df <- data.frame(
      gene = names(ld1_values),
      LD1 = abs(ld1_values)
    )
    # order by absolute LD1 scales
    result_df <- result_df[order(-result_df$LD1),]
    write.csv(result_df, file = paste0("ranked_markers_pos_",stage, ".csv"), row.names = FALSE)
  }
