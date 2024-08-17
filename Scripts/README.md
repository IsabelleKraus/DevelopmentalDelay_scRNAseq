# Mouse Embryonic Stage Prediction using scRNA-seq Data

### Analysis of the Wild-Type (WT) Reference scRNA-seq Data

- **WT_stats.r**  
  Generates general statistics for the WT dataset, including the number of cells, embryos, and the distribution of male and female embryos across the dataset.

- **WT_UMAP.r**  
  Creates Uniform Manifold Approximation and Projection (UMAP) visualizations showing WT single-cell transcriptomes from six mouse embronic stages (E6.5-E9.0) and 55 cell type annotated clusters

- **WT_stacked_plots.r**  
  Produces stacked bar plots to illustrate the percentage of each cell type within individual embryos, highlighting the cellular composition at different stages of development.

### Evaluation and Benchmarking of Machine Learning Methods

The scripts in this section evaluate and benchmark three different methods for predicting embryonic stages based on scRNA-seq data.

- **File Naming Convention**  
  - **\*10.r**: Scripts use 10% of randomly selected WT cells for validation, while the remaining 90% are used for training.
  - **\*3embryos.r**: Scripts use cells from 3 embryos per developmental stage (3x6 embryos in total) for validation, with the remaining WT cells from 54 embryos used for training.

#### PCA on Cell Type Information

- **Devdel_PCA_10.r** and **Devdel_PCA_3embryos.r**  
  Perform Principal Component Analysis (PCA) on cell state information to predict embryonic stages. The PCA is conducted on the median proportions of cell states and their presence/absence across developmental stages. The `prcomp` function is used to compute the principal components, which are then applied to classify validation embryos based on the closest WT stage using Euclidean distance.

#### Label Transfer Using Seurat

- **Devdel_LT_10.r** and **Devdel_LT_3embryos.r**  
  Adapt the Seurat package's Label Transfer method for embryonic stage prediction. The `FindTransferAnchors` function identifies anchors between WT training data and validation data, and the `TransferData` function uses these anchors to transfer stage information.

#### LDA on Stage-Specific Marker Genes

- **Devdel_LDA_10.r** and **Devdel_LDA_3embryos.r**  
  Apply Linear Discriminant Analysis (LDA) to predict embryonic stages based on stage-specific marker genes identified using the Seurat function `FindMarkers`. LDA is performed using the MASS package's `lda` function. The LDA model is then used to classify validation cells and embryos by stage, with predictions based on the posterior probabilities of each stage.

### Marker Gene Analysis

- **Stage_pos_markers.r**  
  Generates ranked and unranked lists of stage-specific marker genes. It also creates an upset plot to visualize how many marker genes are unique to specific stages.

- **GSEA_ORA_stage_pos_markers.r**  
  Performs Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) on stage-specific marker genes using the WebGestaltR package.
For GSEA, marker genes are ranked based on their absolute LD1 values. The script generates bubble plots to visualize the results of both ORA and GSEA, as well as UMAPs to display the expression patterns of highly ranked marker genes.

### Staging of Knockout (KO) Embryos using LDA

- **Devdel_LDA_KO.r**  
  Identifies stage-specific marker genes and conducts a new LDA using the complete set of WT cells as the training reference. This model is then applied to stage the KO cells. The process follows the same LDA approach used for the evaluations, with KO embryos classified based on the single-cell stage classification learned from the WT reference. The script also generates plots that visualize the staging results for the knockout cells and embryos.
