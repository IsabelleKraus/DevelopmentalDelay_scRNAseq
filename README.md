# Quantification of developmental delay in early embryogenesis using scRNAseq

## Abstract 
Embryonic development is a highly coordinated process that can be disrupted in mutants, and quantifying these delays is crucial for understanding the role of regulator genes. This project utilizes single-cell RNA sequencing (scRNA-seq) data to classify embryonic stages and quantify developmental delays in mouse embryos. Three machine learning methods for the stage prediction based on scRNA-seq data were comprehensively evaluated and compared. The first method involved Principal Component Analysis (PCA) on cell state compositions. The second method, Label Transfer, which is typically used to transfer cell type information, was adapted to understand the degree to which stage information can also be learned. Lastly, a novel approach using Linear Discriminant Analysis (LDA) on stage-specific differentially expressed genes (markers) was implemented to predict the staging of single cells. The LDA model demonstrated the highest performance in stage classification and was further applied to knockout (KO) embryos to assess developmental delays. 

## Contents
This repository contains R scripts for:
- The analysis of annotated scRNA-seq data.
- The benchmarking of three methods for the prediction of embryonic stages:
  - PCA on cell type sompositions
  - Label Transfer using Seurat
  - LDA on stage-specific marker gene expression
- The analysis of stage-specific marker genes
- The staging of knockout embryos using the LDA model

## Data availability 
All scRNA-seq Seurat objects, stage-specific marker genes and the LDA model are available at https://nc.molgen.mpg.de/cloud/index.php/s/xyiAcn7PCpb8QTN
