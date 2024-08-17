library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(stringr)

load("WT.Robj")


# Annotations
cluster_order <- c(2, 38, 6, 53, 41, 10, 48, 17, 32, 37, 13, 52, 3, 44, 21, 34, 14, 45, 46, 50, 28, 4, 42, 54, 0, 30, 16, 47, 20, 18, 25, 33, 7, 9, 24, 31, 8, 49, 51, 43, 29, 5, 19, 23, 35, 26, 22, 39, 27, 11, 15, 1, 40, 36, 12)
cluster_names <- c('2_Epiblast', '38_PostProxPrimStreak', '6_EarlyMeso', '53_PGC', '41_BodyWall', '10_ProxExE', '48_DiffExE', '17_DiffExE', '32_Trophoblasts', '37_DefEndo', '13_PrimEndo', '52_YolkSac2', '3_YolkSac1', '44_ParietalEndo', '21_GutEndo', '34_HemogenEndothel', '14_Angioblast', '45_MonocMacrophProg', '46_EMP', '50_HematoEndoProg', '28_PrimBloodProg', '4_Blood_early', '42_PrimFetal_Blood', '54_PrimFetal_Blood', '0_Blood_late', '30_Xmeso_early', '16_Allantois', '47_AmnioticMeso3', '20_AmnioticMeso2', '18_AmnioticMeso1', '25_PresomMeso', '33_EarlyMeso_post', '7_PostLat_IntermMeso', '9_SecHF', '24_PharyArchMeso', '31_PrimHF', '8_Somites', '49_GenitourMeso', '51_Artifact', '43_NodeNotochord', '29_IndEpi_early', '5_IndEpi_late', '19_SurfaceEcto', '23_NeuralRidge_post', '35_NeuralRidge_ant', '26_PrimStreak', '22_NMP', '39_PreneuralPlate_post', '27_PreneuralPlate_ant', '11_Forebrain', '15_Midbrain', '1_NeuralPlate', '40_FloorPlate', '36_MotorNeuron', '12_NeuralCrest')
cluster_colors <- c('#000000', '#474747', '#a0a0a0', '#DB083D', '#722c61', '#fcb7b7', '#f4595c', '#db484a', '#5d0000', '#f6c100', '#f7a800', '#f89e00', '#f98e00', '#fb7100', '#f5f200', '#fad4a6', '#deb88c', '#cba47a', '#bf986f', '#b28b63', '#a57d57', '#98704b', '#8b633f', '#7c5432', '#673f1f', '#2cfeff', '#5bdffe', '#7acafe', '#a9abfd', '#da8afc', '#d767ff', '#a746d0', '#9339bd', '#7927a4', '#62188e', '#3f006c', '#166788', '#099be2', '#492b00', '#1142fc', '#c2ecb3', '#99e2a9', '#4ad094', '#25c88a', '#00bf80', '#00ffda', '#00efca', '#00e1bc', '#00d2ac', '#00c19a', '#00ac83', '#009e75', '#008d62', '#00784b', '#003d0e')
cluster_group <- c('Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Placenta', 'Placenta', 'Placenta', 'Placenta', 'Embryo', 'YolkSac', 'YolkSac', 'YolkSac', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo')
cluster_germlayer <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'Meso', 'XEcto', 'XEcto', 'XEcto', 'XEcto', 'Endo', 'XEndo', 'XEndo', 'XEndo', 'XEndo', 'Endo', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'XEndo', 'Meso', 'Epiblast', 'Epiblast', 'Ecto', 'Ecto', 'Ecto', 'Meso', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto')
# Load the stringr library

# Define the cluster names with IDs and names
cluster_full_names <- c(
  "Epiblast (2)", "Posterior proximal primitive streak (38)", "Early mesoderm (6)", 
  "Primordial germ cells (53)", "Body wall (41)", "Proximal extraembryonic ectoderm (10)", 
  "Differentiated extraembryonic ectoderm 1 (48)", "Differentiated extraembryonic ectoderm 2 (17)", 
  "Trophoblasts (32)", "Definitive endoderm (37)", "Primitive endoderm (13)", "Yolk sac 2 (52)", 
  "Yolk sac 1 (3)", "Parietal endoderm (44)", "Gut endoderm (21)", "Hemogenic endothelium (34)", 
  "Angioblast (14)", "Monocyte/Macrophage progenitor (45)", "Erythro-myeloid progenitor (46)", 
  "Hematopoietic endothelial progenitor (50)", "Primitive blood progenitor (28)", "Early blood (4)", 
  "Primitive fetal blood 1 (54)", "Primitive fetal blood 2 (42)", "Late blood (0)", 
  "Early extraembryonic mesoderm (30)", "Allantois (16)", "Amniotic mesoderm 3 (47)", 
  "Amniotic mesoderm 2 (20)", "Amniotic mesoderm 1 (18)", "Presomitic mesoderm (25)", 
  "Posterior early mesoderm (33)", "Posterior lateral intermediate mesoderm (7)", 
  "Secondary heart field (9)", "Pharyngeal arch mesoderm (24)", "Primary heart field (31)", 
  "Somites (8)", "Genitourinary mesoderm (49)", "Artifact (51)", "Node notochord (43)", 
  "Early induced epiblast (29)", "Late induced epiblast (5)", "Surface ectoderm (19)", 
  "Posterior neural ridge (23)", "Anterior neural ridge (35)", "Primitive streak (26)", 
  "Neuromesodermal progenitor (22)", "Posterior preneural plate (39)", "Anterior preneural plate (27)", 
  "Forebrain (11)", "Midbrain (15)", "Neural plate (1)", "Floor plate (40)", "Motor neuron (36)", 
  "Neural crest (12)"
)

# Extract numbers and names
numbers <- as.integer(sub(".*\\((\\d+)\\)", "\\1", cluster_full_names))
names <- sub("(.*)\\s\\(\\d+\\)", "\\1", cluster_full_names)
# Format with right-aligned numbers and names
cluster_full_names <- paste(str_pad(numbers, width = 2, side = "left", pad = " "), "\t", names, sep=" ")

annotation <- data.frame(State=cluster_order, cluster_colors, cluster_names, cluster_group, cluster_germlayer, cluster_full_names)
rownames(annotation) <- annotation$cluster_names
annotation$State <- as.character(annotation$State)

# colors and cluster_order combined
named_colors <- setNames(cluster_colors, cluster_full_names)
unique_stages <- unique(WT@meta.data$stage)
unique_embryos <- unique(WT@meta.data$embryo)
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')
named_stages <- setNames(e_days, unique_stages)
# color of lineages 
col <- c('#000000', '#a9abfd', '#722c61', '#DB083D',  '#fcb7b7', '#f6c100', '#7c5432', '#00d2ac', '#fb7100'  )
lineages <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'XEcto', 'Endo', 'Blood', 'Ecto', 'XEndo')
germ_colors <- setNames(col,lineages) 

merged_data <- WT@meta.data %>% left_join(annotation, by= c("State"="State"))
WT <- AddMetaData(WT, metadata = merged_data$cluster_germlayer, col.name = "cluster_germlayer")
WT <- AddMetaData(WT, metadata = merged_data$cluster_group, col.name = "cluster_group")
WT <- AddMetaData(WT, metadata = merged_data$cluster_names, col.name = "cluster_names")
WT <- AddMetaData(WT, metadata = merged_data$cluster_colors, col.name = "cluster_colors")
WT <- AddMetaData(WT, metadata = merged_data$cluster_full_names, col.name = "cluster_full_names")


embryo_kernel_cells <- WT@meta.data %>% group_by(stage, embryo, State, sex) %>% summarise(cells_per_state = n()) %>% ungroup() 
total_cells <- embryo_kernel_cells %>% group_by(stage,embryo,sex) %>% summarise(total_cells = sum(cells_per_state)) 
embryo_kernel_cells <- embryo_kernel_cells %>% left_join(total_cells, by = c("stage","embryo", "sex")) %>% mutate(percentage = cells_per_state/total_cells * 100)
embryo_kernel_cells$embryo_num <- gsub("embryo", "", embryo_kernel_cells$embryo)
embryo_kernel_cells <- embryo_kernel_cells %>% left_join(annotation, by= c("State"="State"))
embryo_kernel_cells <- embryo_kernel_cells %>%
        arrange(stage, sex, embryo) %>%
        mutate(order_number = (match(paste(stage, embryo,sex), unique(paste(stage, embryo,sex)))))
embryo_kernel_cells <- embryo_kernel_cells %>% arrange(cluster_germlayer, cluster_full_names)
embryo_kernel_cells$cluster_full_names <- factor(embryo_kernel_cells$cluster_full_names, levels = unique(embryo_kernel_cells$cluster_full_names))
embryo_kernel_cells$sex_label <- ifelse(embryo_kernel_cells$sex == "male", "m", "f")
embryo_sex_label <- embryo_kernel_cells %>%
  group_by(order_number, stage) %>%  # Group by both order_number and stage
  summarize(sex_label = first(sex_label), .groups = 'drop')  # Get the first sex label for each group

# stacked plots with WT reference 
WT_stacked_plot <- ggplot(embryo_kernel_cells, aes(x = factor(order_number), y = percentage, fill = cluster_full_names)) +
        geom_bar(stat = "identity", position = "stack", size=0.5) +
        facet_grid(~stage , scales="free", space="free_x", labeller=as_labeller(named_stages)) +
        theme_minimal() +
        scale_fill_manual(name="Cell State", values= named_colors) +
        labs(x="Embryo", y="Percent of cells")+
        scale_color_manual(name="Cell State", values=named_colors)+
        scale_x_discrete()+
        geom_text(data=embryo_sex_label,aes(x=factor(order_number), label = sex_label, y = 101),vjust = 0.25, inherit.aes = FALSE) +
        theme( axis.text.x=element_text(angle=45),axis.text=element_text(size=16),strip.text = element_text(size = 20), axis.title=element_text(size=16), plot.title = element_text(size = 20 ) )
ggsave("all_WT_stacked.pdf", WT_stacked_plot,units="mm", height=150, width=600 )

# Legend
WT_stacked_plot <- WT_stacked_plot + guides(fill=guide_legend(ncol=5))
l <- get_legend(WT_stacked_plot)
ggsave("cell_stateslegend.pdf", l, width=14)



