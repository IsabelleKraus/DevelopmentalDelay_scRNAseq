# Quantification of developmental delay in early embryogenesis using scRNAseq

## Abstract 
Embryonic development is a highly coordinated process that can be disrupted in mutants, and quantifying these delays is crucial for understanding the role of regulator genes. This project utilizes single-cell RNA sequencing (scRNA-seq) data to classify embryonic stages and quantify developmental delays in mouse embryos. Three machine learning methods for the stage prediction based on scRNA-seq data were comprehensively evaluated and compared. The first method involved Principal Component Analysis (PCA) on cell state compositions. The second method, Label Transfer, which is typically used to transfer cell type information, was adapted to understand the degree to which stage information can also be learned. Lastly, a novel approach using Linear Discriminant Analysis (LDA) on stage-specific differentially expressed genes was implemented to predict the staging of single cells. The LDA model demonstrated the highest performance in stage classification and was further applied to knockout (KO) embryos to assess developmental delays. 

## Contents
This repository contains R scripts used for the analysis of annotated scRNA-seq data, the development of a new method to quantify embryonic staging using an LDA, including a comprehensive benchmark against previously used or widely established methods and its application on knockout embryos.

## Data availability 

