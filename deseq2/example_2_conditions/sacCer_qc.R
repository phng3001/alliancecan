# ############# Initial setup #############

# Set working directory
setwd("/project/def-mouellet/Scripts_MOU/PNP/alliancecan/deseq2/")
getwd()

# Activate environment
renv::activate()

# Load required packages
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

# Set seed
set.seed(42)





# ############# Load data #############

# Gene/transcript count
data <- read.table("example_2_conditions/sacCer_counts_raw.txt", header = T, row.names = 1)

# Metadata
meta <- read.table("example_2_conditions/sacCer_meta.txt", header = T, row.names = 1)

# Convert to factor
meta[] <- lapply(meta, as.factor) # all columns

# Check that sample names match in both files
## Composition check
all(colnames(data) %in% rownames(meta)) # Should be TRUE
## Order check
all(colnames(data) == rownames(meta)) 
## If order check gives FALSE
# data_raw <- data
# data <- data_raw[, match(rownames(meta), colnames(data_raw))]





# ############# QC #############

## ########### Count data modeling ###########
# Poisson: mean == variance; Binomial: mean < variance

# Control condition
mean_counts <- apply(data[, c("WT_C_1", "WT_C_2", "WT_C_3")], 1, mean)
var_counts <- apply(data[, c("WT_C_1", "WT_C_2", "WT_C_3")], 1, var)
df <- data.frame(mean_counts, var_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=var_counts)) +
  geom_line(aes(x=mean_counts, y=mean_counts, color="red"), show.legend = FALSE) +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic()
ggsave("example_2_conditions/control_mean_var_plot.pdf", width = 8, height = 6, dpi = 300) # Extensions: pdf, png, svg etc

# Treatment condition
mean_counts <- apply(data[, c("WT_E_1", "WT_E_2", "WT_E_3")], 1, mean)
var_counts <- apply(data[, c("WT_E_1", "WT_E_2", "WT_E_3")], 1, var)
df <- data.frame(mean_counts, var_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=var_counts)) +
  geom_line(aes(x=mean_counts, y=mean_counts, color="red"), show.legend = FALSE) +
  scale_y_log10() +
  scale_x_log10() +
  theme_classic()
ggsave("example_2_conditions/treatment_mean_var_plot.pdf", width = 8, height = 6, dpi = 300)



## ########### Normalized counts matrix ###########

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~condition
)

# Original count matrix
# counts(dds)

# DESeq2 normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Normalized counts matrix
normalized_counts <- counts(dds, normalized = TRUE)
write.table(normalized_counts, 
            file = "example_2_conditions/normalized_counts.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)



## ########### Normalized counts transformed ###########

### ######### rlog transformation #########
# rlog() works well on small datasets (n < 30)
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)
write.table(rld_mat, 
            file = "example_2_conditions/normalized_counts_rlog.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)

### ######### vst transformation #########
# vst() is faster, recommended for medium to large datasets (n > 30)
vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)
write.table(vsd_mat, 
            file = "example_2_conditions//normalized_counts_vst.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)

### ######### log2(n+1) transformation #########
ntd <- normTransform(dds)
ntd_mat <- assay(ntd)
write.table(ntd_mat, 
            file = "example_2_conditions/normalized_counts_log2.tsv", 
            sep = "\t", quote = FALSE, col.names = NA)



## ########### PCA ###########
# Will use rlog transformation
# vst or log2(n+1) could be used similarly

rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "condition") # Top 500 most variable genes
# plotPCA(rld, intgroup = "condition", ntop = 100) # Top 100 most variable genes

# Save in pdf
pdf("example_2_conditions/PCA.pdf", width = 8, height = 6)
plotPCA(rld, intgroup = "condition")
dev.off()

# Save in png
# png("example_2_conditions/PCA.png", width = 8, height = 6, units = "in", res = 300)
# plotPCA(rld, intgroup = "condition")
# dev.off()

# Save in tiff
# tiff("example_2_conditions/PCA.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
# plotPCA(rld, intgroup = "condition")
# dev.off()

# PCA with ggplot2
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(
  pcaData,
  aes(
    x = PC1, y = PC2,
    color = condition
  )
) +
  geom_point(size = 5) +
  labs(
    x = paste0("PC1 (", percentVar[1], "%)"),
    y = paste0("PC2 (", percentVar[2], "%)")
  ) +
  theme_classic()
ggsave("example_2_conditions/PCA_ggplot.pdf", width = 8, height = 6, dpi = 300)



## ########### Hierarchical Clustering ###########
# Will use rlog transformation
# vst or log2(n+1) could be used similarly

rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)

### ######### Correlation matrix #########

# Pairwise correlation
rld_cor <- cor(rld_mat)

# Default heatmap
pheatmap(rld_cor)

# Customize heatmap
# pheatmap() + pdf() / png() ...  often have problem. May need to re-run sometimes to get output file
heat.colors <- brewer.pal(9, "Blues")
png("example_2_conditions/correlation_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(
  rld_cor, 
  color = heat.colors,
  cluster_cols = TRUE, cluster_rows = TRUE,
  border_color=NA, fontsize = 10, height=20
)
dev.off()



### ######### Count matrix #########
normalized_counts <- counts(dds, normalized = TRUE)

#### ####### Genes highly expressed #######
# Top 20 most expressed genes
top20_expression_genes <- order(rowMeans(normalized_counts), decreasing = TRUE)[1:20]

# Heatmap
png("example_2_conditions/top_20_expression_genes_heatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(
  rld_mat[top20_expression_genes, ],
  cluster_cols = TRUE, cluster_rows = TRUE,
  fontsize_row = 10,
  annotation_col = meta
)
dev.off()

#### ####### Genes highly varied #######
# Top 20 genes with most variance
top20_variance_genes <- order(rowVars(normalized_counts), decreasing = TRUE)[1:20]

# Heatmap
png("example_2_conditions/top_20_variance_genes_heatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(
  rld_mat[top20_variance_genes, ],
  cluster_cols = TRUE, cluster_rows = TRUE,
  fontsize_row = 10,
  annotation_col = meta
)
dev.off()


