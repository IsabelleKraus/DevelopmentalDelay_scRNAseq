library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse) 
library(viridis)
library(WebGestaltR)
library(rbioapi)
library(biomaRt)
library(clusterProfiler)
library(org.Mm.eg.db)
library(data.table)
library(ComplexUpset)

load("WT.Robj")
stages <- c("WT_65", "WT_70", "WT_75", "WT_80", "WT_85", "WT_90")
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, stages)

# Define a set of distinct colors for each stage (same as in PCA)
stage_colors <- rev(viridis(6, option = "C", begin = 0.0, end = 1))
color_mapping <- setNames(stage_colors, e_days)
color_palette <- setNames(stage_colors, c("6.5","7","7.5","8","8.5","9"))


# Import individual stage-specific marker gene lists from CSV, loading only first column
all_markers <- setNames(lapply(stages, function(set_name) {
  # first column contains the genes 
  df <- read.csv(paste0("ranked_markers_pos_", set_name, ".csv"),
                 header = TRUE, colClasses = c("character", "numeric"))
  df
}), stages)

for(s in stages){
  df <- all_markers[[s]]
  write.table(df, file = paste0("ranked_markers_pos_", s, ".rnk"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

####################################################################################################################################
# Gene Set Enrichment Analysis (GSEA) using LD1 ranked genes
####################################################################################################################################

enriched <- list()

# Perform molecular function analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("GSEA_MF_min3")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  # submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "GSEA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Molecular_Function_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGeneFile =  paste0("/project/PRCcomp/data/MarkerGenes/ranked_markers_pos_", stage_name, ".rnk"),
    interestGeneType = "genesymbol",
    minNum = 3,
    sigMethod = "top",
    outputDirectory = paste0("GSEA_MF_min3"),
    projectName = paste0(stage_name,"_GSEA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Perform cellular component analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("GSEA_CC_min3")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  # submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "GSEA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Cellular_Component_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGeneFile =  paste0("ranked_markers_pos_", stage_name, ".rnk"),
    interestGeneType = "genesymbol",
    minNum = 3,
    sigMethod = "top",
    outputDirectory = paste0("GSEA_CC_min3"),
    projectName = paste0(stage_name,"_GSEA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Perform biological process analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("GSEA_BP_min3")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  # submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "GSEA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Biological_Process_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGeneFile =  paste0("ranked_markers_pos_", stage_name, ".rnk"),
    interestGeneType = "genesymbol",
    minNum = 3,
    sigMethod = "top",
    outputDirectory = paste0("GSEA_BP_min3"),
    projectName = paste0(stage_name,"_GSEA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Generate bubble plot for GSEA results (top 5 for each) 
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

  file <- read.table(paste0("GSEA_MF_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t", nrows=5,header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR, size = file$size, type = "MF")
  file <- read.table(paste0("GSEA_BP_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t",  nrows=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR,size = file$size,  type = "BP")
  file <- read.table(paste0("GSEA_CC_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t",  nrows=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR, size = file$size, type = "CC")

  combined <- bind_rows(result_bp, result_cc, result_mf)
  print(min(combined$FDR))

  bubble <- ggplot(combined, aes(x = enrichmentScore, y = reorder(description, enrichmentScore), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red") +  # Adjust colors as needed
    scale_size(name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "GO term", title=stage_name) +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text = element_text(size = 16, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 16) ,
    panel.border = element_blank(),
    axis.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16) )
  ggsave(paste0(stage_name, "_GSEA_bubble_top5.pdf"), bubble, height=10, width=14)
}

# Generate bubble plot for GSEA results 
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

  file <- read.table(paste0("GSEA_MF_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t", header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR, size = file$size, type = "MF")
  file <- read.table(paste0("GSEA_BP_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t", header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR,size = file$size,  type = "BP")
  file <- read.table(paste0("GSEA_CC_min3/Project_", stage_name, "_GSEA_top/enrichment_results_", stage_name,"_GSEA_top.txt"),sep="\t", header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentScore = file$enrichmentScore, pValue = file$pValue, FDR = file$FDR, size = file$size, type = "CC")

  combined <- bind_rows(result_bp, result_cc, result_mf)

  bubble <- ggplot(combined, aes(x = enrichmentScore, y = reorder(description, enrichmentScore), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red") +  # Adjust colors as needed
    scale_size(range = c(1,13),,breaks = c(min(combined$size), max(combined$size)), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "GO term", title=stage_name) +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text.x = element_text(size = 24),  # Increase x axis tick text size
    axis.text.y = element_text(size = 24),  # Increase y axis tick text size
    strip.text = element_text(size = 24) )
  ggsave(paste0(stage_name, "_GSEA_bubble.pdf"), bubble, height=20,  width=14)
}

#####################################################################################################################################################################
# Over Represenatation Analysis (ORA)
#####################################################################################################################################################################

# Perform molecular function analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("ORA_MF")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  data <- fread( paste0("ranked_markers_pos_", stage_name, ".rnk"))
  genes <- as.list(data[[1]])
  # submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Molecular_Function_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGene = genes ,
    referenceSet = "genome",
    interestGeneType = "genesymbol",
    sigMethod = "top",
    outputDirectory = paste0("ORA_MF"),
    projectName = paste0(stage_name,"_ORA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Perform biological process analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("ORA_BP")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  data <- fread( paste0("ranked_markers_pos_", stage_name, ".rnk"))
  genes <- as.list(data[[1]])
  ## submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Biological_Process_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGene = genes ,
    referenceSet = "genome",
    interestGeneType = "genesymbol",
    sigMethod = "top",
    outputDirectory = paste0("ORA_BP"),
    projectName = paste0(stage_name,"_ORA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Perform cellular component analysis for each marker set
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

    output_dir <- paste0("ORA_CC")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

  data <- fread( paste0("ranked_markers_pos_", stage_name, ".rnk"))
  genes <- as.list(data[[1]])
  # submit enrichment request
  gsea_result <- WebGestaltR(
    enrichMethod = "ORA",
    organism = "mmusculus",
    enrichDatabase = "geneontology_Cellular_Component_noRedundant",
    enrichDatabaseType = "genesymbol",
    interestGene = genes ,
    referenceSet = "genome",
    interestGeneType = "genesymbol",
    sigMethod = "top",
    outputDirectory = paste0("ORA_CC"),
    projectName = paste0(stage_name,"_ORA_top"),
    saveRawGseaResult = TRUE
  )
} 

# Generate bubble plots (all results) 
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")

  combined <- bind_rows(result_bp, result_cc, result_mf)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red") +  
    #scale_size_continuous(breaks = seq(min(combined$size), max(combined$size), length.out=5), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "GO term", title=stage_name) +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 12),
    legend.title = element_text(size=12)) 
  ggsave(paste0(stage_name, "_ORA_bubble.pdf"), bubble, height=14, width=10)
}

stage_name <- "WT_65"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red") + 
    scale_size_continuous(range = c(3, 10),breaks = seq(min(combined$size), max(combined$size), by=7), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "GO term", title="") +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text = element_text(size = 12, color="black"),  
    strip.text = element_text(size = 12),
    legend.title = element_text(size=12)) 
  ggsave(paste0(stage_name, "_ORA_bubble.pdf"), bubble, height=10, width=10)

stage_name <- "WT_85"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name= "p Val") +
    scale_size_continuous(range = c(1, 8),breaks = seq(min(combined$size), max(combined$size), by=2), name="Gene count") +  
    labs(x = "Fold enrichment", y = "GO term", title="") +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text = element_text(size = 12, color="black"), 
    strip.text = element_text(size = 12),
    legend.title = element_text(size=12)) 
  ggsave(paste0(stage_name, "_ORA_bubble.pdf"), bubble, height=10, width=10)


# Buble plots with FDR filtered GO terms
stage_name <- "WT_85"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_mf, result_cc)
  combined85 <- combined %>% filter(FDR < 0.05) 

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name="p-value") +  # Adjust colors as needed
    scale_size_continuous(range = c(3,10),breaks=seq(min(combined$size), max(combined$size), by=2), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "", title="") +
    theme_classic()+
    facet_grid(rows=vars(type), scales="free_y")+
    theme(legend.position = "right",
    axis.text = element_text(size = 16, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 16) ,
    panel.border = element_blank(),
    axis.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16))
  ggsave(paste0(stage_name, "_ORA_bubble_top5.pdf"), bubble, height=8, width=14)

stage_name <- "WT_75"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)
  combined75 <- combined %>% filter(FDR < 0.05)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name="p-value") +  # Adjust colors as needed
    scale_size_continuous(range = c(2,9),breaks=seq(min(combined$size), max(combined$size), by=2), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "", title="") +
    theme_classic()+
    facet_grid(rows=vars(type), scales="free_y")+
    theme(legend.position = "right",
    axis.text = element_text(size = 16, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 16) ,
    panel.border = element_blank(),
    axis.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16))
  ggsave(paste0(stage_name, "_ORA_bubble_top5.pdf"), bubble, height=8, width=14)

stage_name <- "WT_65"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)
  combined65 <- combined %>% filter(FDR < 0.05)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name="p-value") +  # Adjust colors as needed
    scale_size_continuous(range = c(3,15),breaks=seq(min(combined$size), max(combined$size), by=8), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "", title="") +
    theme_classic()+
    facet_grid(rows=vars(type), scales="free_y")+
    theme(legend.position = "right",
    axis.text = element_text(size = 16, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 16) ,
    panel.border = element_blank(),
    axis.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16))
  ggsave(paste0(stage_name, "_ORA_bubble_top5.pdf"), bubble, height=8, width=14)

stage_name <- "WT_70"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)
  combined70 <- combined %>% filter(FDR < 0.05)

stage_name <- "WT_80"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)
  combined80 <- combined %>% filter(FDR < 0.05)

stage_name <- "WT_90"
  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t",nrow=5, header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", nrow=5, header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)
  combined90 <- combined %>% filter(FDR < 0.05)

combined65$stage <- "E6.5"
combined75$stage <- "E7.5"
combined85$stage <- "E8.5"
combined90$stage <- "E9.0"
combined80$stage <- "E8.0"
combined70$stage <- "E7.0"

# Bubble plot showing for stages E6.5, E7.5 and E8.5 FDR filtered GO terms 
combined <- bind_rows(combined65,combined75,combined85)
combined <- combined %>% filter(type=="BP")
  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name="p-value") +  # Adjust colors as needed
    scale_size_continuous(range = c(3,15),breaks=seq(min(combined$size), max(combined$size), by=8), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "", title="") +
    theme_classic()+
    facet_grid(rows=vars(stage), scales="free_y")+
    theme(legend.position = "right",
    axis.text = element_text(size = 16, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 16) ,
    panel.border = element_blank(),
    axis.title = element_text(size=16),
    legend.text = element_text(size=16),
    legend.title = element_text(size=16))
  ggsave(paste0("65_75_85", "_ORA_bubble_top5.pdf"), bubble, height=8, width=14)


# Bubble plot showing for each stage FDR filtered GO terms (biological processes)
combined <- bind_rows(combined65,combined70,combined75,combined80,combined85,combined90)
combined <- combined %>% filter(type=="BP")
  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = pValue)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue", high = "red", name="p-value") +  # Adjust colors as needed
    scale_size_continuous(range = c(3,15),breaks=seq(min(combined$size), max(combined$size), by=8), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "Biological process", title="") +
    theme_classic()+
    facet_grid(rows=vars(stage), scales="free_y", space="free")+
    theme(legend.position = "right",
    panel.border = element_rect(colour = "black", fill=NA),
    axis.text = element_text(size = 20, color="black"),  # Increase x axis tick text size
    strip.text = element_text(size = 20), 
    axis.title = element_text(size=20),
    legend.text = element_text(size=20),
    legend.title = element_text(size=20))
  ggsave(paste0("all_BP", "_ORA_bubble_top5.pdf"), bubble, height=14, width=14)


# Bubble plot for each stage colored by FDR values
for (stage_name in names(all_markers)) {
  cat("Processing:", stage_name, "\n")

  file <- read.table(paste0("ORA_MF/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_mf <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "MF")
  file <- read.table(paste0("ORA_BP/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_bp <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR,size = file$overlap,  type = "BP")
  file <- read.table(paste0("ORA_CC/Project_", stage_name, "_ORA_top/enrichment_results_", stage_name,"_ORA_top.txt"),sep="\t", header=TRUE)
  result_cc <- data.frame(description=file$description, enrichmentRatio = file$enrichmentRatio, pValue = file$pValue, FDR = file$FDR, size = file$overlap, type = "CC")
  combined <- bind_rows(result_bp, result_cc, result_mf)

  bubble <- ggplot(combined, aes(x = enrichmentRatio, y = reorder(description, enrichmentRatio), size = size, color = FDR)) +
    geom_point(alpha=0.7, show.legend = TRUE) +
    scale_color_gradient(low = "blue",high = "grey") +  # Adjust colors as needed
    scale_size(range = c(1,13), name="Gene count") +   # Adjust bubble sizes and legend title
    labs(x = "Fold enrichment", y = "GO term", title=stage_name) +
    facet_grid(rows=vars(type), scales="free_y")+
    theme_minimal() +
    theme(legend.position = "right",
    axis.text.x = element_text(size = 12),  # Increase x axis tick text size
    axis.text.y = element_text(size = 12),  # Increase y axis tick text size
    strip.text = element_text(size = 14) )
  ggsave(paste0(stage_name, "_ORA_bubble_FDR.pdf"), bubble, height=14, width=10)
}

########################################################################################################################################################################
# UMAPs of wild-type reference showing expression values of highly LD1-ranked marker genes  
########################################################################################################################################################################

plot <- FeaturePlot(WT,  features= "Dppa5a" , order=TRUE, pt.size=1)+NoAxes()
ggsave("umap_LDA_Dppa5a.pdf", plot)

plot <- FeaturePlot(WT,  features= "Hbb-y" , order=TRUE, pt.size=1)+NoAxes()
ggsave("umap_LDA_Hbby.pdf", plot)

plot <- FeaturePlot(WT,  features= "Elf5" , order=TRUE, pt.size=1)+NoAxes()
ggsave("umap_LDA_Elf5.pdf",plot)

plot <- FeaturePlot(WT,  features= "Pou5f1" , order=TRUE, pt.size=1)+NoAxes()
ggsave("umap_LDA_Pou5f1.pdf",plot)

plot <- FeaturePlot(WT,  features= c("Dppa5a","Elf5") ,order=TRUE, pt.size=1)
ggsave("umap_LDA_Dppa5a_Elf5.pdf",plot,width=14)

########################################################################################################################################################################
# Overlap of marker genes for different wild-type training sets 
########################################################################################################################################################################

em_marker <- read.csv(file="pos_marker_genes_3embryos.csv", header=TRUE)$Gene
ten_marker <- read.csv(file="pos_marker_genes_10.csv", header=TRUE)$Gene
WT_marker <- read.csv(file="pos_marker_genes_WT.csv", header=TRUE)$Gene

length(intersect(em_marker, WT_marker))
length(intersect(em_marker, WT_marker))/length(WT_marker) # 96.7% overlapping 

length(intersect(em_marker, ten_marker))
length(intersect(em_marker, ten_marker))/length(ten_marker) # 96.3% overlapping 

length(intersect(WT_marker, ten_marker))
length(intersect(WT_marker, ten_marker))/length(WT_marker) # 99.3% overlapping 

########################################################################################################################################################################
# Overlap of marker genes in different enriched GO terms of E6.5 (Upset plot) 
########################################################################################################################################################################

mf <- data.frame(read.table("ORA_MF/Project_WT_65_ORA_top/enrichment_results_WT_65_ORA_top.txt",sep="\t", header=TRUE))
bp <- data.frame(read.table("ORA_BP/Project_WT_65_ORA_top/enrichment_results_WT_65_ORA_top.txt",sep="\t", header=TRUE))
cc <- data.frame(read.table("ORA_CC/Project_WT_65_ORA_top/enrichment_results_WT_65_ORA_top.txt",sep="\t", header=TRUE))

create_gene_df <- function(data) {
data %>%
  tidyr::separate_rows(userId, sep = ";") %>%
  dplyr::rename(Gene = userId) %>%
  dplyr::select(Gene, description)
}

# data frames for each category
bp_df <- create_gene_df(bp)
mf_df <- create_gene_df(mf)
cc_df <- create_gene_df(cc)
bp_df$source <- "bp"
mf_df$source <- "mf"
cc_df$source <- "cc"
combined_df <- bind_rows(bp_df, mf_df, cc_df)
gene_data_wide <- combined_df %>%
  mutate(value = 1) %>%
  spread(key = description, value = value, fill = 0)
gene_data_binary <- gene_data_wide %>%
  dplyr::select(-Gene, -source)

# Upset Plot showing intersection sizes >2, without singletons 
upset <- upset(gene_data_binary, names(gene_data_binary),name="GO term intersection",min_size=2, min_degree=2, height_ratio=1.2,themes=upset_default_themes(text=element_text(size=24)))
ggsave("Upset_GO_terms_65.pdf", upset, width=19, height=14)








