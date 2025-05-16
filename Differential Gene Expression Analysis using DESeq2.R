# Differential Gene Expression Analysis using DESeq2

# Install Required Packages (if not already installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  BiocManager::install("DESeq2")
}
if (!requireNamespace("RUVSeq", quietly = TRUE)) {
  BiocManager::install("RUVSeq")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  BiocManager::install("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  BiocManager::install("RColorBrewer")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# Load Packages
library(DESeq2)
library(RUVSeq)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Import Gene Counts
COUNTS <- read.csv(file="C:/Users/Lenovo/Downloads/GSE227516_counts.csv", header=TRUE, row.names=1)
head(COUNTS)

## -----------------------------------------------------------------------------
# Import Sample Information
META <- read.csv(file="C:/Users/Lenovo/Downloads/sample_information.csv", header=TRUE)
head(META)

## -----------------------------------------------------------------------------
# Reorder columns
COUNTS <- COUNTS[, c("P1", paste0("P", 2:9), "P10")]
head(COUNTS)

## -----------------------------------------------------------------------------
# Rounding off and convert to matrix
COUNTS <- round(COUNTS)
COUNTS <- as.matrix(COUNTS)

## -----------------------------------------------------------------------------
# List unique sample conditions
unique(META$condition)

## ----warning=FALSE, message=FALSE---------------------------------------------
# Creating DESeq2 Dataset
dds <- DESeqDataSetFromMatrix(countData = COUNTS, 
                              colData=META, 
                              design=~condition)

## -----------------------------------------------------------------------------
dim(dds)

## -----------------------------------------------------------------------------
# Filtering low count genes
threshold <- 10
dds <- dds[ rowMeans(counts(dds)) >= threshold,]

## ----warning=FALSE, message=FALSE---------------------------------------------
# DESeq2 Analysis
prdds <- DESeq(dds)
prdds


## -----------------------------------------------------------------------------
# Normalization
norm_counts <- counts(prdds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)
head(norm_counts)

## -----------------------------------------------------------------------------
# Transformation
mks <- estimateSizeFactors(dds)
rld <- rlogTransformation(prdds, blind = FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

## ----scatter_plots, dev='png', fig.show='hide'--------------------------------
# Scatter Plots Comparison
par(mfrow=c(1, 3))
lims <- c(-2, 20)
plot(log2(counts(mks, normalized=TRUE)[,1:2] + 1),pch=16, cex=0.3, main="log2(x + 1)", xlim=lims, ylim=lims)
plot(assay(rld)[,1:2], pch=16, cex=0.3, main="R log", xlim=lims, ylim=lims)
plot(assay(vsd)[,1:2], pch=16, cex=0.3, main="VST", xlim=lims, ylim=lims)


## ----hist_plots, dev='png', fig.show='hide'-----------------------------------
# Histograms Comparison
par(mfrow=c(1, 3))
hist(counts(mks))
hist(assay(rld))
hist(assay(vsd))


## ----s2s_heatmap_plot, dev='png', fig.show='hide'-----------------------------
# Sample-to-sample distances
sample_dist <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dist)
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sample_dist_matrix,
         clustering_distance_rows=sample_dist,
         clustering_distance_cols=sample_dist,
         col=colors,)

## ----pca_plot, dev='png', fig.show='hide', message=FALSE, warning=FALSE-------
# PCA Plot
pca_data <- plotPCA(rld, intgroup = c("condition"), returnData = TRUE)
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 2, aes(color = condition)) +
  geom_text_repel(aes(label = rownames(pca_data)), nudge_x = 0, nudge_y = 0) +
  xlab(paste0("PC1: ", round(attr(pca_data, "percentVar")[1], 2) * 100, "% variance")) +
  ylab(paste0("PC2: ", round(attr(pca_data, "percentVar")[2], 2) * 100, "% variance"))


## ----dispersion_plot, dev='png', fig.show='hide'------------------------------
# Dispersion Plot
plotDispEsts(prdds, main = "Dispersion plot", 
             genecol="gray20", fitcol="red", 
             finalcol="dodgerblue3" 
) 

## -----------------------------------------------------------------------------
# DESeq2 Result
res05 <- results(prdds, alpha = 0.05)
res05 <- na.omit(res05)

## -----------------------------------------------------------------------------
# Order by adjusted p-value
res05ordered <- res05[order(res05$padj),]
head(as.data.frame(res05ordered))


## ----ma_plot, dev='png', fig.show='hide'--------------------------------------
# MA Plot
DESeq2::plotMA(
  res05, 
  main="Sedentary vs Exercise, alpha=0.05", 
  ylim=c(-5,10),
  cex=0.5, 
  colNonSig=adjustcolor("gray20", alpha.f=0.5), 
  colSig=adjustcolor("dodgerblue3", alpha.f=0.5) 
)
abline(h = 1, col = '#ff0000' , lwd = 1)
abline(h = -1, col= '#ff0000', lwd = 1)


## ----volcano_plot, dev='png', fig.show='hide'---------------------------------
# Volcano Plot
res05$gene_status <- ifelse(
  res05$padj < 0.05, 
  ifelse(
    res05$log2FoldChange > 1, 
    "Up-Regulated",
    ifelse(
      res05$log2FoldChange < -1, 
      "Down-Regulated", 
      "Non-significant"
    )
  ), 
  "Non-significant"
)

ggplot(
  res05, 
  aes(x = log2FoldChange, y = -log10(padj), color = factor(gene_status))
) +
  geom_point(size = 1, alpha = 0.7) +
  scale_color_manual(values = brewer.pal(3, "Set1")) +
  theme_minimal() +
  ggtitle("Volcano Plot of Differentially Expressed Genes") +
  xlab("log2 Fold Change") +
  ylab("-log10(Adjusted p-value)") +
  theme(legend.title = element_blank()) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed")

## -----------------------------------------------------------------------------
sig_genes <- as.data.frame(res05[res05$padj < 0.05 & abs(res05$log2FoldChange) > 1, ])
head(sig_genes)


## -----------------------------------------------------------------------------
up_genes <- subset(sig_genes, log2FoldChange > 0)
head(up_genes)

## -----------------------------------------------------------------------------
down_genes <- subset(sig_genes,log2FoldChange < 0)
head(down_genes)

## ----upreg_heatmap, dev='png', fig.show='hide'--------------------------------
top_up <- head(up_genes[order(up_genes$log2FoldChange, decreasing = TRUE), ], 10)
top_up_exp <- assay(rld)[rownames(top_up), ]
pheatmap(top_up_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "Top Up Regulated Genes Heatmap")


## ----downreg_heatmap, dev='png', fig.show='hide'------------------------------
top_down <- head(down_genes[order(down_genes$log2FoldChange, decreasing = TRUE), ], 10)
top_down_exp <- assay(rld)[rownames(top_down), ]
pheatmap(top_down_exp,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "row",
         show_colnames = TRUE,
         col=brewer.pal(name="RdBu", n=11),
         main = "Top Down Regulated Genes Heatmap")


## -----------------------------------------------------------------------------
n_total <- nrow(COUNTS)
n_de_genes <- nrow(sig_genes)
n_up_genes <- nrow(up_genes)
n_down_genes <- nrow(down_genes)

## -----------------------------------------------------------------------------
# Combine significant genes with their categories
sig_genes$gene_status <- ifelse(sig_genes$log2FoldChange > 0, "Up-Regulated", "Down-Regulated")
sig_genes$gene_id <- rownames(sig_genes)

# Write the data frame to a file
if (!file.exists("result")) {
  dir.create("result")
}

write.table(sig_genes, file = "result/significant_DE_genes.csv", sep = ",", row.names = FALSE)


## -----------------------------------------------------------------------------
sessionInfo()
