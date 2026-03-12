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
library(ggrepel)
library(clusterProfiler)
library(org.Sc.sgd.db)

# Set seed
set.seed(42)





# ############# Load data #############

## ########### Count data ###########

data <- read.table("example_2_conditions/input/sacCer_counts_raw.txt", 
                   header = T, row.names = 1)
str(data)

# If count data is not integer, round to integer
# data <- data %>%
#   mutate(across(where(is.numeric), ~ as.integer(round(.))))



## ########### Col data ###########

meta <- read.table("example_2_conditions/input/sacCer_meta.txt", 
                   header = T, row.names = 1)

# Convert to factor
meta[] <- lapply(meta, as.factor) # all columns
str(meta)

## ########### Check sample names match ###########

## Composition check
all(colnames(data) %in% rownames(meta)) # Should be TRUE

## Order check
all(colnames(data) == rownames(meta)) 

## If order check gives FALSE
# data_raw <- data
# data <- data_raw[, match(rownames(meta), colnames(data_raw))]





# ############# QC #############

## ########### Count distribution ###########

# Data long format
data_long <- data %>%
  gather(
    colnames(data),
    key = "sample_name",
    value = "raw_counts"
  )

# Count distribution
ggplot(
  data_long,
  aes(x = raw_counts)
) +
  geom_histogram(stat = "bin", bins = 100) +
  facet_wrap(~ sample_name) +
  xlab("Raw expression counts") +
  ylab("Number of genes") +
  theme_bw()
# ggsave("example_2_conditions/output/count_distribution_plot.pdf", 
#        width = 8, height = 6, dpi = 300) # Extensions: pdf, png, svg etc



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
# ggsave("example_2_conditions/output/control_mean_var_plot.pdf", 
#        width = 8, height = 6, dpi = 300)



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
# ggsave("example_2_conditions/output/treatment_mean_var_plot.pdf", 
#        width = 8, height = 6, dpi = 300)



## ########### Normalized counts matrix ###########

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~condition
)

# Original count matrix
head(counts(dds))

# DESeq2 normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

# Normalized counts matrix
normalized_counts <- counts(dds, normalized = TRUE)
# write.table(normalized_counts, 
#             file = "example_2_conditions/output/normalized_counts.tsv", 
#             sep = "\t", quote = FALSE, col.names = NA)



## ########### Normalized counts transformed ###########

### ######### rlog transformation #########
# rlog() works well on small datasets (n < 30)
rld <- rlog(dds, blind = TRUE)
rld_mat <- assay(rld)
# write.table(rld_mat, 
#             file = "example_2_conditions/output/normalized_counts_rlog.tsv", 
#             sep = "\t", quote = FALSE, col.names = NA)

### ######### vst transformation #########
# vst() is faster, recommended for medium to large datasets (n > 30)
vsd <- vst(dds, blind = TRUE) # n genes >= 1000
# vsd <- varianceStabilizingTransformation(dds, blind = TRUE) # n genes < 1000
vsd_mat <- assay(vsd)
# write.table(vsd_mat, 
#             file = "example_2_conditions/output/normalized_counts_vst.tsv", 
#             sep = "\t", quote = FALSE, col.names = NA)

### ######### log2(n+1) transformation #########
ntd <- normTransform(dds)
ntd_mat <- assay(ntd)
# write.table(ntd_mat, 
#             file = "example_2_conditions/output/normalized_counts_log2.tsv", 
#             sep = "\t", quote = FALSE, col.names = NA)



## ########### PCA ###########
# Will use rlog transformation
# vst or log2(n+1) could be used similarly

# rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "condition") # Top 500 most variable genes
# plotPCA(rld, intgroup = "condition", ntop = 100) # Top 100 most variable genes

# Save in pdf
# pdf("example_2_conditions/output/PCA.pdf", width = 8, height = 6)
plotPCA(rld, intgroup = "condition")
# dev.off()

# Save in png
# png("example_2_conditions/output/PCA.png", width = 8, height = 6, units = "in", res = 300)
# plotPCA(rld, intgroup = "condition")
# dev.off()

# Save in tiff
# tiff("example_2_conditions/output/PCA.tiff", width = 8, height = 6, units = "in", res = 300, compression = "lzw")
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
# ggsave("example_2_conditions/output/PCA_ggplot.pdf", width = 8, height = 6, dpi = 300)



## ########### Hierarchical Clustering ###########
# Will use rlog transformation
# vst or log2(n+1) could be used similarly

rld_mat <- read.table("example_2_conditions/output/normalized_counts_rlog.tsv", 
                      header = T, row.names = 1)
head(rld_mat)

### ######### Correlation matrix #########

# Pairwise correlation
rld_cor <- cor(rld_mat)

# Default heatmap
pheatmap(rld_cor)

# Customize heatmap
# pheatmap() + pdf() / png() on Rstudio ...  often have problem. May need to re-run sometimes to get the right output file
heat.colors <- brewer.pal(9, "Blues")
# png("example_2_conditions/output/correlation_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(
  rld_cor, 
  color = heat.colors,
  cluster_cols = TRUE, cluster_rows = TRUE,
  border_color=NA, fontsize = 10, height=20
)
# dev.off()



### ######### Count matrix #########

normalized_counts <- read.table("example_2_conditions/output/normalized_counts.tsv",
                                header = T, row.names = 1)
head(normalized_counts)

#### ####### Genes highly expressed #######

# Top 20 most expressed genes
top20_expression_genes <- order(rowMeans(normalized_counts), decreasing = TRUE)[1:20]
normalized_counts[top20_expression_genes,] %>%
  row.names()

# Heatmap
# png("example_2_conditions/output/top_20_expression_genes_heatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(
  rld_mat[top20_expression_genes, ],
  cluster_cols = TRUE, cluster_rows = TRUE,
  fontsize_row = 10,
  annotation_col = meta
)
# dev.off()



#### ####### Genes highly varied #######

# Top 20 genes with most variance
top20_variance_genes <- order(rowVars(as.matrix(normalized_counts)), decreasing = TRUE)[1:20]
normalized_counts[top20_variance_genes,] %>%
  row.names()

# Heatmap
# png("example_2_conditions/output/top_20_variance_genes_heatmap.png", width = 8, height = 6, units = "in", res = 300)
pheatmap(
  rld_mat[top20_variance_genes, ],
  cluster_cols = TRUE, cluster_rows = TRUE,
  fontsize_row = 10,
  annotation_col = meta
)
# dev.off()





# ############# Differential expression analysis #############

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~condition
)

# Original count matrix
head(counts(dds))

## ########### Pre-filter (optional) ###########
# https://support.bioconductor.org/p/65256/

# Option 1: Filter genes where there are less than 10 counts
# idx <- rowSums(counts(dds)) >= 10
# dds <- dds[idx,]

# Option 2: Filter genes where there are less than 3 samples having at least 5 counts
# idx <- rowSums(counts(dds) >= 5) >= 3
# dds <- dds[idx,]

# Option 3: Filter genes where there are less than 3 samples with normalized counts greater than or equal to 5
# dds <- estimateSizeFactors(dds)
# idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
# dds <- dds[idx,]

# Option 4: Only keep genes that have a total count of at least 10 in any condition
# cond <- colData(dds)$condition
# idx <- rep(FALSE, nrow(dds))
# for (c in levels(cond)) {
#   idx <- idx | rowSums(counts(dds)[, cond == c, drop=FALSE]) >= 10
# }
# dds <- dds[idx,]

# Option 5: Only keep genes that have at least 10 count in every sample of at least 1 condition
# cond <- colData(dds)$condition
# idx <- rep(FALSE, nrow(dds))
# for (c in levels(cond)) {
#   mat <- counts(dds)[, cond == c, drop = FALSE]
#   idx <- idx | rowSums(mat >= 10) == ncol(mat)
# }
# dds <- dds[idx,]



## ########### Run DESeq2 ###########

# Set reference level
dds$condition <- relevel(dds$condition, ref = "C")

# Run DE analysis
dds <- DESeq(dds)



### ######### Dispersion estimates #########
# Dispersion (α): Var = μ + α*μ^2
# Dispersion ~ 1/μ, ~ Var

# png("example_2_conditions/output/dispersion_estimates.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds)
# dev.off()



### ######### Result exploration #########

# Available coefficients
resultsNames(dds)

#### ####### Unshrunken LFC #######

# alpha: adjusted p-value cutoff (FDR)
res <- results(
  dds, 
  name = "condition_E_vs_C", 
  alpha = 0.05
  )

# Summarize results
summary(res)

# MA plot
plotMA(res)

# Convert to dataframe, sort by padj
res_df <- res %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  arrange(padj)

# Save to file
# write.table(res_df, 
#             file = "example_2_conditions/output/treatment_vs_control_unshrunken_result.tsv", 
#             sep = "\t", quote = FALSE, row.names = FALSE)



#### ####### Shrunken LFC #######
# Shrinkage of LFC estimates toward zero when the information for a gene is low (low counts, high dispersion) 
# This does not change the total number of genes identified as significantly DE

# apeglm shrinkage: only for use with 'coef'
resLFC <- lfcShrink(
  dds, 
  coef = "condition_E_vs_C", 
  res = res, 
  type = "apeglm"
  )

# Summarize results
summary(resLFC)

# MA plot
plotMA(resLFC)

# Convert to dataframe, sort by padj
resLFC_df <- resLFC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  arrange(padj)

# Save to file
# write.table(resLFC_df, 
#             file = "example_2_conditions/output/treatment_vs_control_result.tsv", 
#             sep = "\t", quote = FALSE, row.names = FALSE)

# Column information
mcols(resLFC, use.names = TRUE)



#### ####### DE genes #######
# Use shrunken LFC

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- log2(2)

# Significant DE genes
resLFC_df_sig <- resLFC_df %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Save to file
write.table(resLFC_df_sig, 
            file = "example_2_conditions/output/treatment_vs_control_result_DE.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)



### ######### Result visualisation #########

# Normalized count data
normalized_counts <- read.table("example_2_conditions/output/normalized_counts.tsv",
                                header = T, row.names = 1)
normalized_counts_df <- normalized_counts %>%
  rownames_to_column(var = "gene")
head(normalized_counts_df)

# Metadata
meta_df <- meta %>%
  rownames_to_column("sample_name")
meta_df

#### ####### Count plot #######

##### ##### A single gene #####

plotCounts(dds, gene = "YDL248W", intgroup = "condition")

normalized_counts_df %>%
  filter(gene == "YDL248W")



##### ##### Multiple genes #####

# Example: Top 20 DE genes by padj
top20_sig_genes <- resLFC_df %>%
  arrange(padj) %>%
  pull(gene) %>%
  head(n=20)
top20_sig_genes

# Subset normalized counts for top 20 genes
top20_sig_norm <- normalized_counts_df %>%
  filter(gene %in% top20_sig_genes)

# Convert to long format
top20_sig_norm <- top20_sig_norm %>%
  gather(
    colnames(top20_sig_norm)[-1], 
    key = "sample_name", 
    value = "normalized_counts"
  )

# Merge with metadata
top20_sig_norm <- inner_join(meta_df, top20_sig_norm, by = "sample_name")

# Boxplots
ggplot(
  top20_sig_norm,
  aes(x = gene, y = normalized_counts, color = condition)
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  scale_color_manual(
    labels = c("Control", "Treatment"),
    values = c("gray", "blue")
  ) +
  labs(
    x = NULL, 
    y = "log10 normalized counts", 
    title = "Top 20 Significant DE Genes"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
# ggsave("example_2_conditions/output/treatment_vs_control_top_20_DE_boxplot.pdf", 
#        width = 8, height = 6, dpi = 300)

# Scatterplots
ggplot(
  top20_sig_norm, 
  aes(x = condition, y = normalized_counts, color = condition)
) +
  geom_jitter(height = 0, width = 0.15) +
  scale_y_continuous(trans = 'log10') +
  scale_color_manual(
    labels = c("C", "E"),
    values = c("gray", "blue")
  ) +
  ylab("log10 normalized expression level") +
  xlab("condition") +
  ggtitle("Top 20 Significant DE Genes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ gene) +
  theme_bw()
# ggsave("example_2_conditions/output/treatment_vs_control_top_20_DE_scatterplot.pdf", 
#        width = 8, height = 6, dpi = 300)



#### ####### Heatmap #######

# Subset normalized counts for significant genes
sig_norm <- normalized_counts_df %>%
  filter(gene %in% resLFC_df_sig$gene) %>%
  column_to_rownames("gene")

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Heatmap
# png("example_2_conditions/output/treatment_vs_control_DE_heatmap.png", width = 8, height = 8, units = "in", res = 300)
pheatmap(
  sig_norm,
  color = heat_colors,
  cluster_rows = TRUE, cluster_cols = TRUE,
  show_rownames = FALSE,
  annotation_col = meta,
  fontsize = 10,
  border_color = NA,
  scale = "row",
  fontsize_row = 10,
  height = 20
)
# dev.off()
# Notes: Normalized gene counts are scaled by z-score standardization



#### ####### Volcano plot #######

# resLFC_df <- read.table("example_2_conditions/output/treatment_vs_control_result.tsv", 
#                         sep = "\t", header = T)

# Set thresholds
# padj.cutoff <- 0.05
# lfc.cutoff <- log2(2)

# Specify up/down/ns expressed genes
resLFC_df_annot <- resLFC_df %>%
  mutate(
    expression = case_when(
      padj < padj.cutoff & log2FoldChange > lfc.cutoff ~ "up",
      padj < padj.cutoff & log2FoldChange < -lfc.cutoff ~ "down",
      TRUE ~ "ns"
    )
  )

# Data range
## x axis
range(resLFC_df_annot$log2FoldChange, na.rm = TRUE)
## y axis
max(-log10(range(resLFC_df_annot$padj, na.rm = TRUE)))



#### ##### Simple Volcano plot #####

resLFC_df_annot %>%
  ggplot(
    aes(x = log2FoldChange, y = -log10(padj), color = expression, alpha = 0.5)
  ) +
  geom_point(size = 3) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(-lfc.cutoff, lfc.cutoff),
    linetype = "dashed"
  ) +
  scale_color_manual(
    breaks = c("up", "down", "ns"),
    values = c("tomato", "steelblue", "grey")
  ) +
  guides(alpha = "none") +
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 200)) +
  labs(
    title = "Gene expression changes between Treatment and Control",
    x = "log2(FC)", 
    y = "-log10(padj)", 
    color = "Expression\nchange"
  ) +
  theme_classic()



#### ##### Volcano plot with top genes labeling #####

##### ##### Example 1: Top 10 |log2FC| #####

resLFC_df_annot_label <- resLFC_df_annot %>%
  filter(expression != "ns") %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  arrange(desc(abs_log2FoldChange)) %>%
  mutate(label = "")
num_keep <- 10
resLFC_df_annot_label$label[1:num_keep] <- resLFC_df_annot_label$gene[1:num_keep]

resLFC_df_annot %>%
  ggplot(
    aes(x = log2FoldChange, y = -log10(padj), color = expression, alpha = 0.5)
  ) +
  geom_point(size = 3) +
  geom_text_repel(
    data = resLFC_df_annot_label,
    aes(label = label),
    max.overlaps = Inf,
    alpha = 1
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = c(-lfc.cutoff, lfc.cutoff),
    linetype = "dashed"
  ) +
  scale_color_manual(
    breaks = c("up", "down", "ns"),
    values = c("tomato", "steelblue", "grey")
  ) +
  guides(alpha = "none") +
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 200)) +
  labs(
    title = "Gene expression changes between Treatment and Control",
    x = "log2(FC)", y = "-log10(padj)", color = "Expression\nchange"
  ) +
  theme_classic()



##### ##### Example 2: Selection based on |log2FC| and padj #####

resLFC_df_annot %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  filter(-log10(padj) > 50 & abs(log2FoldChange) > 5) %>%
  dim()

resLFC_df_annot_label <- resLFC_df_annot %>%
  mutate(abs_log2FoldChange = abs(log2FoldChange)) %>%
  filter(-log10(padj) > 50 & abs(log2FoldChange) > 5) %>%
  mutate(label = gene)

ggplot(
  data = resLFC_df_annot,
  aes(x = log2FoldChange, 
      y = -log10(padj))
) +
  geom_point(
    aes(colour = expression),
    alpha = 0.6,
    shape = 19,
    size = 2
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) + 
  geom_vline(
    xintercept = c(-lfc.cutoff, lfc.cutoff),
    linetype = "dashed"
  ) +
  geom_label_repel(
    data = resLFC_df_annot_label,
    aes(label = label),
    force = 2
  ) +
  scale_color_manual(
    breaks = c("up", "down", "ns"),
    values = c("darkred", "darkblue", "grey")
  ) +
  scale_x_continuous(
    breaks = c(seq(-8, 8, 2)),
    limits = c(-8, 8)
  ) +
  labs(
    title = "Gene expression changes between Treatment and Control",
    x = "log2(FC)", 
    y = "-log10(padj)", 
    color = "Expression\nchange"
  ) +
  theme_light()
# ggsave("example_2_conditions/output/treatment_vs_control_DE_volcano.pdf", 
#        width = 8, height = 6, dpi = 300)





# ############# Functional analysis #############

## ########### Gene lists ###########

# resLFC_df <- read.table("example_2_conditions/output/treatment_vs_control_result.tsv", 
#                         sep = "\t", header = T)

# Background dataset (all genes tested for significance)
all_genes <- as.character(resLFC_df$gene)
length(all_genes)

# Significant DE genes
sig_genes <- as.character(
  resLFC_df %>%
    filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
    pull(gene)
)
length(sig_genes)

# Up & down regulated genes
up_genes <- as.character(
  resLFC_df %>%
    filter(padj < padj.cutoff & log2FoldChange > lfc.cutoff) %>%
    pull(gene)
)

down_genes <- as.character(
  resLFC_df %>%
    filter(padj < padj.cutoff & log2FoldChange < -lfc.cutoff) %>%
    pull(gene)
)



## ########### Over-representation analysis (ORA) ###########

### ######### Up regulated genes #########

ora_bp_up <- enrichGO(
  gene = up_genes,        # Gene list of interest
  universe = all_genes,   # Background genes
  OrgDb = org.Sc.sgd.db,  # OrgDb for S. cerevisiae
  keyType = "ORF",        # Key type of gene identifiers in gene list
  ont = "BP",             # Ontology: BP (Biological Process)
  pAdjustMethod = "BH",   # Multiple testing correction
  pvalueCutoff  = 0.05,   # P-value cutoff
  qvalueCutoff  = 0.2     # Q-value cutoff
)

# Output results to a dataframe
ora_bp_up_df <- data.frame(ora_bp_up)
range(ora_bp_up_df$p.adjust, na.rm = TRUE)

# Save to file
# write.table(ora_bp_up_df, 
#             file = "example_2_conditions/output/treatment_vs_control_result_DE_up_ORA_BP.tsv", 
#             sep = "\t", quote = FALSE, row.names = FALSE)

# Top n terms by p.adjust
# Display gene ratio: Nb of genes related to GO term / total number of significant genes
dotplot(
  ora_bp_up, 
  showCategory = 10, 
  title="GO Biological Processes Enrichment in Up-regulated genes"
  )

# Top n terms by p.adjust
# Display gene count
# png("example_2_conditions/output/treatment_vs_control_result_DE_up_ORA_BP_top10.png", 
#     width = 8, height = 6, units = "in", res = 300)
barplot(
  ora_bp_up, 
  showCategory = 10, 
  title="GO Biological Processes Enrichment in Up-regulated genes"
  )
# dev.off()

ora_bp_up_df %>%
  arrange(p.adjust) %>%
  head(n=10) %>%
  pull(Description)



### ######### Down regulated genes #########

ora_bp_down <- enrichGO(
  gene = down_genes,      # Gene list of interest
  universe = all_genes,   # Background genes
  OrgDb = org.Sc.sgd.db,  # OrgDb for S. cerevisiae
  keyType = "ORF",        # Key type of gene identifiers
  ont = "BP",             # Ontology: BP (Biological Process)
  pAdjustMethod = "BH",   # Multiple testing correction
  pvalueCutoff  = 0.05,   # P-value cutoff
  qvalueCutoff  = 0.2     # Q-value cutoff
)

# Output results to a dataframe
ora_bp_down_df <- data.frame(ora_bp_down)
range(ora_bp_down_df$p.adjust, na.rm = TRUE)

# Save to file
# write.table(ora_bp_down_df, 
#             file = "example_2_conditions/output/treatment_vs_control_result_DE_down_ORA_BP.tsv", 
#             sep = "\t", quote = FALSE, row.names = FALSE)

# Top n terms by p.adjust
# Display gene ratio: Nb of genes related to GO term / total number of significant genes
dotplot(
  ora_bp_down, 
  showCategory = 10, 
  title="GO Biological Processes Enrichment in Down-regulated genes"
  )

# Top n terms by p.adjust
# Display gene count
# png("example_2_conditions/output/treatment_vs_control_result_DE_down_ORA_BP_top10.png", 
#     width = 8, height = 6, units = "in", res = 300)
barplot(
  ora_bp_down, 
  showCategory = 10, 
  title="GO Biological Processes Enrichment in Down-regulated genes"
  )
# dev.off()

ora_bp_down_df %>%
  arrange(p.adjust) %>%
  head(n=10) %>%
  pull(Description)



### ######### Parameters to consider / play with #########  
# DE thresholds: padj, log2FoldChange
# Gene list: DE genes, Up, Down
# Ontology: "BP", "MF", "CC" or "ALL"
# Ontology pAdjustMethod
# Ontology pvalueCutoff, qvalueCutoff

