library(ggplot2)
library(dplyr)
library(Seurat)
library(tidyr)
library(RColorBrewer)

load("WT.Robj")
# To get the object names:
ls()
# Should show that it's a Seurat Object: 
class(WT)
# Seurat Objects have multi slot structure:
names(WT)
# Gene Expression Data: 
head(WT@assays$RNA@counts[, 1:5])
# Meta Data
head(WT@meta.data)

# Dimension of matrix (first num shoudld be the amount of single cells)
dim(WT@meta.data)

unique(WT@meta.data$orig.ident)


# Distribution of metadata:

# How many males and females?
table(WT$sex)

# How many different embryos in total?
length(unique(WT@meta.data$embryo))

# How many empryos per time point?
embryo_per_stage <- table(WT@meta.data$stage, WT@meta.data$embryo)
count_unique_embryo <- apply(embryo_per_stage, 1, function(row) sum(row>0))
count_unique_embryo

# Amount of (fe)males
sex_per_stage <- table(WT@meta.data$stage, WT@meta.data$sex)
sex_per_stage

pdf("WT_sex_distr.pdf")
barplot(table(WT@meta.data$sex),ylab="Frequency")
dev.off()

# Sequenced cells:
mean_cells_per_embryo <- apply(embryo_per_stage, 1, mean)
mean_cells_per_embryo

total_cells <- apply(embryo_per_stage, 1, sum)
total_cells

# To find out how many male and female embryos there are per stage:

female_embryos_per_stage <- WT@meta.data %>% group_by(stage, embryo) %>% filter(sex=="female")%>%
  summarise() %>% group_by(stage) %>%
  summarise(num_females = n()) %>%
  arrange()
female_embryos_per_stage

male_embryos_per_stage <- WT@meta.data %>% group_by(stage, embryo) %>% filter(sex=="male")%>%
  summarise() %>% group_by(stage) %>%
  summarise(num_males = n()) %>%
  arrange()
male_embryos_per_stage

# To find out how many percent of each embryo were represented in cells

embryo_percentage_per_stage <- WT@meta.data %>%
  group_by(stage, embryo) %>%
  summarise(embryo_count = n()) %>%                       
  mutate(total_count = sum(embryo_count),                 
  percentage = (embryo_count / total_count) * 100) %>% 
  select(-embryo_count, -total_count) %>%                 
  arrange(stage, embryo)
embryo_percentage_per_stage

max_embryo_percentage <- embryo_percentage_per_stage %>% group_by(stage) %>% summarise(max = max(percentage)) %>% arrange(stage)
min_embryo_percentage <- embryo_percentage_per_stage %>% group_by(stage) %>% summarise(min = min(percentage)) %>% arrange(stage)
mean_embryo_percentage <- embryo_percentage_per_stage %>% group_by(stage) %>% summarise(mean = mean(percentage)) %>% arrange(stage)
mean_embryo_percentage

# Visualize QC metrics as a violin plot
pdf("QC_Vln_plots")
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()


grouped_meta_data <- WT@meta.data %>% group_by(stage)

# percent of cells (for every kernel) over single embryos 
# single plot for each time point  

# Calculate Percentage of cells:

stages <- unique(WT@meta.data$stage)
embryo_kernel_cells_combined <- tibble()
for (stage in stages) {
    # Filter data for current stage 
    data_subset <- WT@meta.data[WT@meta.data$stage == stage, ]
    #print(head(data_subset))
    # Calculate percentage of cells per embryo and Kernel

    embryo_kernel_cells <- data_subset %>% group_by(embryo, Kernel, sex) %>% summarise(cells_per_kernel = n()) %>% ungroup() 
    total_cells <- embryo_kernel_cells %>% group_by(embryo, sex) %>% summarise(total_cells = sum(cells_per_kernel)) 
    embryo_kernel_cells <- embryo_kernel_cells %>% left_join(total_cells, by = c("embryo", "sex")) %>% mutate(percentage = cells_per_kernel/total_cells * 100)
    embryo_kernel_cells$stage = stage

    embryo_kernel_cells$embryo_num <- as.integer(gsub("embryo", "", embryo_kernel_cells$embryo))
    embryo_kernel_cells <- embryo_kernel_cells %>% arrange(embryo_num)

    embryo_kernel_cells_combined <- rbind(embryo_kernel_cells_combined, embryo_kernel_cells)
}


n_colors <- length(unique(embryo_kernel_cells_combined$Kernel))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
my_colors <- sample(col_vector, n_colors)
named_colors <- setNames(my_colors, levels(embryo_kernel_cells_combined$Kernel))


# Stacked Plot

pdf("WT_stacked_Kernel_new_plots", width=14, height=6)
ggplot(embryo_kernel_cells_combined, aes(x = as.factor(embryo_num), y = percentage, fill = Kernel)) +
    geom_bar(stat = "identity", position = "stack") +
    facet_wrap(~ stage) +
    theme_minimal() +
    scale_fill_manual(name="Kernel", values= named_colors) +
    labs(title = "Percentage of Cells per Embryo and Kernel by Stage",x = "Embryo", y = "Percentage") +
    scale_x_discrete()+
    theme(legend.position = "bottom")
dev.off()



# Summarize the data to get the associated sex .
summarized_data <- embryo_kernel_cells_combined %>%
  group_by(embryo_num) %>%
  summarise(sex = first(sex))

# Counts for Genes and cells
counts_per_cell <- Matrix::colSums(WT@assays$RNA@counts)
counts_per_gene <- Matrix::rowSums(WT@assays$RNA@counts)
genes_per_cell <- Matrix::colSums(WT@assays$RNA@counts>0) # only counting single matrixcells which are true (>0)
cells_per_gene <- Matrix::rowSums(WT@assays$RNA@counts>0)

# Plots: log10 scaled x-Axis, counts + 1 to avoid log(0)
pdf("WT_counts_per_cell.pdf")
ylimits <- c(0,21000)
hist(log10(counts_per_cell+1), ylim=ylimits, main='counts per cell',col='darkolivegreen3')
dev.off()

pdf("WT_genes_per_cell.pdf")
hist(log10(genes_per_cell+1), ylim=ylimits, main='genes per cell',col='darkolivegreen3')
dev.off()


pdf("WT_counts_per_gene.pdf")
hist(log10(counts_per_gene+1), main='counts per gene',col='darkolivegreen3')
dev.off()

pdf("WT_counts_vs_genes_per_cell.pdf")
plot(log10(counts_per_cell+1), log10(genes_per_cell+1), main= "counts vs genes per cell",col='darkolivegreen3')
dev.off()


# UMAP Plot:

umap_data <- as.data.frame(WT[["umap"]]@cell.embeddings)
head(umap_data)


# used umap data which is already in WT@umap:
# and coloured umap data based on Kernels
# so Kernels are assigned to UMAP based on barcodes

umap_data$Kernel <- WT@meta.data[row.names(umap_data),"Kernel"]

pdf("/project/PRCcomp/images/WT_UMAP_plot.pdf")
ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = factor(Kernel))) +
  geom_point(alpha = 0.6, size = 1) +
  theme_minimal() +
  labs(title = "UMAP Plot grouped by Kernel", x = "UMAP_1", y = "UMAP_2") +
  scale_color_manual(name = "Kernel", values = named_colors)
dev.off()

# how to colour UMAP based on same lineage? 
# how to assign lineage to kernels?


umap_data$stage <- WT@meta.data[row.names(umap_data), "stage"]


# Percentage of mitochondrial RNA per stage 
mean_mt_per_stage <- WT@meta.data %>% group_by(stage) %>% summarize(mean_mt_percent = mean(percent.mt)) %>% arrange()
# Maximum 
max(WT@meta.data$percent.mt)


# Scatter Plot number of genes against number of transcripts 
pdf("WT_feature_vs_count_scatter.pdf")
FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
  + geom_abline(intercept = 0, slope = 1, linetype="dashed") 
  + ylim(0,40000) 
  + xlim(0,40000)
dev.off()

#  Scatter Plot percentage of mtRNA against number of transcripts
pdf("WT_mt_vs_count_scatter.pdf")
FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "percent.mt") 
  + geom_abline(intercept = 0, slope = 1, linetype="dashed") 
dev.off()

# PCA Plot
pdf("WT_PCA_plot.pdf")
DimPlot(WT, reduction = "pca")
dev.off()

