library(ggpubr)
library(dplyr)
library("grid") 
library("gridExtra") 
library(Seurat)
library(ggplot2)
library(tidyr)
library(farver)
library(reshape2)

load("WT.Robj")

cluster_order <- c(2, 38, 6, 53, 41, 10, 48, 17, 32, 37, 13, 52, 3, 44, 21, 34, 14, 45, 46, 50, 28, 4, 42, 54, 0, 30, 16, 47, 20, 18, 25, 33, 7, 9, 24, 31, 8, 49, 51, 43, 29, 5, 19, 23, 35, 26, 22, 39, 27, 11, 15, 1, 40, 36, 12)
cluster_names <- c('2_Epiblast', '38_PostProxPrimStreak', '6_EarlyMeso', '53_PGC', '41_BodyWall', '10_ProxExE', '48_DiffExE', '17_DiffExE', '32_Trophoblasts', '37_DefEndo', '13_PrimEndo', '52_YolkSac2', '3_YolkSac1', '44_ParietalEndo', '21_GutEndo', '34_HemogenEndothel', '14_Angioblast', '45_MonocMacrophProg', '46_EMP', '50_HematoEndoProg', '28_PrimBloodProg', '4_Blood_early', '42_PrimFetal_Blood', '54_PrimFetal_Blood', '0_Blood_late', '30_Xmeso_early', '16_Allantois', '47_AmnioticMeso3', '20_AmnioticMeso2', '18_AmnioticMeso1', '25_PresomMeso', '33_EarlyMeso_post', '7_PostLat_IntermMeso', '9_SecHF', '24_PharyArchMeso', '31_PrimHF', '8_Somites', '49_GenitourMeso', '51_Artifact', '43_NodeNotochord', '29_IndEpi_early', '5_IndEpi_late', '19_SurfaceEcto', '23_NeuralRidge_post', '35_NeuralRidge_ant', '26_PrimStreak', '22_NMP', '39_PreneuralPlate_post', '27_PreneuralPlate_ant', '11_Forebrain', '15_Midbrain', '1_NeuralPlate', '40_FloorPlate', '36_MotorNeuron', '12_NeuralCrest')
cluster_colors <- c('#000000', '#474747', '#a0a0a0', '#DB083D', '#722c61', '#fcb7b7', '#f4595c', '#db484a', '#5d0000', '#f6c100', '#f7a800', '#f89e00', '#f98e00', '#fb7100', '#f5f200', '#fad4a6', '#deb88c', '#cba47a', '#bf986f', '#b28b63', '#a57d57', '#98704b', '#8b633f', '#7c5432', '#673f1f', '#2cfeff', '#5bdffe', '#7acafe', '#a9abfd', '#da8afc', '#d767ff', '#a746d0', '#9339bd', '#7927a4', '#62188e', '#3f006c', '#166788', '#099be2', '#492b00', '#1142fc', '#c2ecb3', '#99e2a9', '#4ad094', '#25c88a', '#00bf80', '#00ffda', '#00efca', '#00e1bc', '#00d2ac', '#00c19a', '#00ac83', '#009e75', '#008d62', '#00784b', '#003d0e')
cluster_group <- c('Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Placenta', 'Placenta', 'Placenta', 'Placenta', 'Embryo', 'YolkSac', 'YolkSac', 'YolkSac', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'YolkSac', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo', 'Embryo')
cluster_germlayer <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'Meso', 'XEcto', 'XEcto', 'XEcto', 'XEcto', 'Endo', 'XEndo', 'XEndo', 'XEndo', 'XEndo', 'Endo', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'Blood', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'XMeso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'Meso', 'XEndo', 'Meso', 'Epiblast', 'Epiblast', 'Ecto', 'Ecto', 'Ecto', 'Meso', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto', 'Ecto')

annotation <- data.frame(State=cluster_order, cluster_colors, cluster_names, cluster_group, cluster_germlayer)
rownames(annotation) <- annotation$cluster_names
annotation$State <- as.character(annotation$State)
merged_data <- WT@meta.data %>% left_join(annotation, by= c("State"="State"))

# colors and cluster_order combined
named_colors <- setNames(cluster_colors, cluster_order)
unique_stages <- unique(WT@meta.data$stage)

# colors of germlayers 
col <- c('#000000', '#a9abfd', '#722c61', '#DB083D',  '#fcb7b7', '#f6c100', '#7c5432', '#deb88c', '#fb7100'  )
germ <- c('Epiblast', 'XMeso', 'Meso', 'PGC', 'XEcto', 'Endo', 'Blood', 'Ecto', 'XEndo')
germ_colors <- setNames(col,germ) 
e_days <- c('E6.5', 'E7.0', 'E7.5', 'E8.0', 'E8.5', 'E9.0')


# UMAP Plot

umap_data <- as.data.frame(WT[["umap"]]@cell.embeddings)
head(umap_data)

umap_data$State <- WT@meta.data[row.names(umap_data),"State"]
umap_data$stage <- WT@meta.data[rownames(umap_data), "stage"]
umap_data <- umap_data %>% left_join(annotation, by=c("State"="State"))

plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = cluster_names) )+
  geom_point(size = 1) +
  theme_minimal() +
  #labs(title = "UMAP") +
  theme_void() +
  scale_color_manual(name = "State", values = named_colors)+
  guides(color="none")
ggsave("WT_UMAP_plot_size1.pdf", plot)

# Legend: 
named_colors <- setNames(cluster_colors, cluster_names)

names_plot <- ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = factor(cluster_names))) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "UMAP") +
  theme_void() +
  scale_color_manual(name = "State", values = named_colors)
  #guides(color="none")
legend <- get_legend(names_plot)
ggsave(filename = "Legend_States.pdf", plot = legend)

cluster_full_names <- c("Epiblast", "Posterior Proximal Primitive Streak", "Early Mesoderm", "Primordial Germ Cells", "Body Wall", "Proximal Extraembryonic Ectoderm", "Differentiated Extraembryonic Ectoderm 1","Differentiated Extraembryonic Ectoderm 2", "Trophoblasts", "Definitive Endoderm", "Primitive Endoderm", "Yolk Sac 2", "Yolk Sac 1", "Parietal Endoderm", "Gut Endoderm", "Hemogenic Endothelium", "Angioblast", "Monocyte/Macrophage Progenitor", "Erythro-Myeloid Progenitor", "Hematopoietic Endothelial Progenitor", "Primitive Blood Progenitor", "Early Blood Cells", "Primitive Fetal Blood 1","Primitive Fetal Blood 2", "Late Blood", "Early XMesoderm", "Allantois", "Amniotic Mesoderm 3", "Amniotic Mesoderm 2", "Amniotic Mesoderm 1", "Presomitic Mesoderm", "Posterior Early Mesoderm", "Posterior Lateral Intermediate Mesoderm", "Secondary Heart Field", "Pharyngeal Arch Mesoderm", "Primary Heart Field", "Somites", "Genitourinary Mesoderm", "Artifact", "Node Notochord", "Early Induced Epiblast", "Late Induced Epiblast", "Surface Ectoderm", "Posterior Neural Ridge", "Anterior Neural Ridge", "Primitive Streak", "Neuromesodermal Progenitor", "Posterior Preneural Plate", "Anterior Preneural Plate", "Forebrain", "Midbrain", "Neural Plate", "Floor Plate", "Motor Neuron", "Neural Crest")
df <- data.frame(Name= cluster_full_names, Color=cluster_colors)
named_cluster_names <- setNames(cluster_colors,cluster_full_names)

df$NewColumn <- 1
names_plot <- ggplot(df, aes(x = 1, y = 1, fill = Name)) +
  geom_tile() +
  #scale_color_manual(name="Cell States", values= named_cluster_names) +
  theme_minimal() +
  scale_fill_manual(name="Cell State", values= named_cluster_names) +
  scale_color_manual(values=named_cluster_names)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
legend <- get_legend(names_plot)
ggsave(filename = "Legend_States_b.pdf", plot = legend, width=12)


# UMAP Plots per Stage

named_colors <- setNames(cluster_colors, cluster_order)

# Create the UMAP plots by stage
pdf("WT_UMAP_stages_wrapped_plot.pdf", width=18, height=3)
facet_plot <- ggplot(umap_data, aes(x= UMAP_1, y=UMAP_2, color = factor(State))) +
  geom_point(size=0.5) +
  theme_minimal() +
  labs(title= "UMAP per Stage", x= "UMAP_1", y="UMAP_2") +
  scale_color_manual(name="Kernel", values = named_colors) +
  facet_wrap(~ stage, ncol = length(unique_stages))+
  theme_void() +
  guides(color="none")
  # if you want a legend within, delete guides() and uncomment theme()
  #theme(legend.position = "bottom")
facet_plot
dev.off()

# Legend: 
named_colors <- setNames(cluster_colors, cluster_names)

# UMAP Plot per Stage with all stages underlying in grey:
umap_data_stage <- dplyr::select(umap_data, -stage) 

stage.labs = e_days
names(stage.labs) <- unique_stages

#pdf("/project/PRCcomp/images/WT_UMAP_stages_wrapped_plot_test.pdf", width=18, height=3)
facet_plot <- ggplot(umap_data, aes(x= UMAP_1, y=UMAP_2)) +
  geom_point(data=umap_data_stage,aes(x=UMAP_1, y=UMAP_2), size=1, colour="grey82", alpha=0.5) +
  theme_minimal() +
  scale_color_manual(name="State", values = named_colors) +
  facet_wrap(~ stage, nrow = length(unique_stages), labeller= labeller(stage=stage.labs), strip.position="left")+
  geom_point(aes(color=cluster_names), size=1) + 
  theme(
        panel.border = element_rect(colour = "black",fill=NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 28, angle = 0)) +
  guides(color = "none")
#facet_plot
#dev.off()
ggsave("WT_UMAP_stages_vs_all_plot_b.pdf", facet_plot, height=6*length(unique_stages), width=6.5)


#pdf("/project/PRCcomp/images/WT_UMAP_stages_wrapped_plot_test.pdf", width=18, height=3)
facet_plot <- ggplot(umap_data, aes(x= UMAP_1, y=UMAP_2)) +
  geom_point(data=umap_data_stage,aes(x=UMAP_1, y=UMAP_2), size=1, colour="grey82", alpha=0.5) +
  theme_minimal() +
  scale_color_manual(name="State", values = named_colors) +
  facet_wrap(~ stage, ncol = length(unique_stages), labeller= labeller(stage=stage.labs), strip.position="top")+
  geom_point(aes(color=cluster_names), size=1) + 
  theme(
        panel.border = element_rect(colour = "black",fill=NA, linewidth = 1),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 48, angle = 0)) +
  guides(color = "none")
#facet_plot
#dev.off()
ggsave("WT_UMAP_stages_vs_all_plot_horizontal.pdf", facet_plot, width=6*length(unique_stages), height=6.5)

