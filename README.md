# Differential-Gene-Expression-Analysis-using-DESeq2
This repository contains R scripts and guidance for performing Differential Gene Expression (DGE) analysis using the [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) package. It is designed to help researchers identify significantly differentially expressed genes from RNA-Seq data.


## 📌 Features

- Import and preprocess count data
- Construct DESeq2 dataset
- Perform normalization and transformation
- Run differential expression analysis
- Visualize results: MA plot, volcano plot, PCA, heatmaps
- Export results to CSV

- ## 🧬 Requirements

- R (>= 4.0)
- Bioconductor packages:
  - `DESeq2`
  - `pheatmap`
  - `ggplot2`
  - `EnhancedVolcano`
  - `apeglm` (optional for shrinkage)
- CSV/TSV file with raw counts (genes × samples)
- Metadata file (sample information, including condition/groups)


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "pheatmap", "EnhancedVolcano", "apeglm"))
install.packages("ggplot2")


Prepare your data:

counts_matrix.csv: Rows are genes, columns are sample counts.

metadata.csv: Contains sample info including condition/group column.


Output
deseq2_results.csv: Table of DE genes with log2FC, p-values, adjusted p-values

Plots:

MA plot

Volcano plot

Heatmap of top genes

PCA of samples


![Volcano Plot of Differentially Expressed Genes](https://github.com/user-attachments/assets/249b68bc-2ea2-4349-9ff8-85a8c6d237b8)



![Top Up Regulated Genes Heatmap](https://github.com/user-attachments/assets/d36aa5fd-d457-4a0f-bd5f-3c710cd736cb)



![Top Down Regulated Genes Heatmap](https://github.com/user-attachments/assets/66e9b068-78ec-406e-b9ed-f7a37f4e7685)






Follow like share and subscribe https://www.youtube.com/@Bioinformatics_Made_Easy

