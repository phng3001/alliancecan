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





# ############# Differential expression analysis #############

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~condition
)

# Original count matrix
# counts(dds)
# str(counts(dds))

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



## ########### Check the model fit ###########

# Dispersion (α): Var = μ + α*μ^2
# Dispersion ~ 1/μ, ~ Var
png("example_2_conditions/dispersion_estimates.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds)
dev.off()



## ########### Explore the results ###########

# Available coefficients
resultsNames(dds)

### ######### Unshrunken LFC #########

# alpha: adjusted p-value cutoff (FDR)
res <- results(dds, 
               name = "condition_E_vs_C", 
               alpha = 0.05)

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
write.table(res_df, 
            file = "example_2_conditions/treatment_vs_control_unshrunken_result.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)



### ######### Shrunken LFC #########
# Shrinkage of LFC estimates toward zero when the information for a gene is low (low counts, high dispersion) 
# This does not change the total number of genes identified as significantly DE

# apeglm shrinkage: only for use with 'coef'
resLFC <- lfcShrink(dds, 
                    coef = "condition_E_vs_C", 
                    res = res, 
                    type = "apeglm")

# Summarize results
summary(resLFC)

# MA plot
plotMA(resLFC)

# Column information
mcols(resLFC, use.names = TRUE)

# Convert to dataframe, sort by padj
resLFC_df <- resLFC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  arrange(padj)

# Save to file
write.table(resLFC_df, 
            file = "example_2_conditions/treatment_vs_control_result.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)



### ######### DE genes #########
# Use shrunken LFC

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- log2(2)

# Significant DE genes
resLFC_df_sig <- resLFC_df %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Save to file
write.table(resLFC_df_sig, 
            file = "example_2_conditions/treatment_vs_control_result_DE.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)



## ########### Visualisation ###########

# Add row names to metadata
metadata <- meta %>%
  rownames_to_column("sample_name")

# normalized_counts <- counts(dds, normalized = TRUE)
normalized_counts <- read.table("example_2_conditions/normalized_counts.tsv", 
                                sep = "\t", header = T, row.names = 1)

# Convert row names to a column named "gene"
normalized_counts <- normalized_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene")



### ######### Count plot #########

#### ####### A single gene #######

plotCounts(dds, gene = "YDL248W", intgroup = "condition")

normalized_counts %>%
  filter(gene == "YDL248W")



#### ####### Top DE genes #######

# Top 20 DE genes by padj
top20_sig_genes <- resLFC_df %>%
  arrange(padj) %>%
  pull(gene) %>%
  head(n=20)

# Subset normalized counts for top 20 genes
top20_sig_norm <- normalized_counts %>%
  filter(gene %in% top20_sig_genes)

# Convert to long format
top20_sig_norm_long <- top20_sig_norm %>%
  gather(
    colnames(top20_sig_norm)[2:7], 
    key = "sample_name", 
    value = "normalized_counts"
    )

# Merge with metadata
top20_sig_norm_long_with_meta <- inner_join(metadata, top20_sig_norm_long, by = "sample_name")

##### ##### Boxplot #####
ggplot(
  top20_sig_norm_long_with_meta,
  aes(x = gene, y = normalized_counts, color = condition)
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  labs(
    x = NULL, y = "log10 normalized counts", title = "Top 20 Significant DE Genes"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

ggplot(
  top20_sig_norm_long_with_meta,
  aes(x = gene, y = normalized_counts, color = condition)
) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  scale_y_log10() +
  scale_color_manual(
    labels = c("Control", "Treatment"),
    values = c("gray", "blue")
  ) +
  labs(
    x = NULL, y = "log10 normalized counts", title = "Top 20 Significant DE Genes"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.title = element_text(hjust = 0.5)
  )
ggsave("example_2_conditions/treatment_vs_control_top_20_DE_boxplot.pdf", width = 8, height = 6, dpi = 300)

# plotCounts(dds, gene = "YPL223C", intgroup = "condition")



##### ##### Scatterplot #####

ggplot(
  top20_sig_norm_long_with_meta, 
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
ggsave("example_2_conditions/treatment_vs_control_top_20_DE_scatterplot.pdf", width = 8, height = 6, dpi = 300)



### ######### Heatmap #########

# Subset normalized counts for significant genes
sig_norm <- normalized_counts %>%
  filter(gene %in% resLFC_df_sig$gene) %>%
  column_to_rownames("gene")

# Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")

# Heatmap
# Sometimes need to re-run to get the output file
png("example_2_conditions/treatment_vs_control_DE_heatmap.png", width = 8, height = 8, units = "in", res = 300)
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
dev.off()
# Notes: Normalized gene counts are scaled by z-score standardization



### ######### Volcano plot #########

resLFC_df <- read.table("example_2_conditions/treatment_vs_control_result.tsv", 
                        sep = "\t", header = T)

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- log2(2)

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

#### ####### Simple Volcano plot #######
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



#### ####### Volcano plot with top genes labeling #######

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
ggsave("example_2_conditions/treatment_vs_control_DE_volcano.pdf", width = 8, height = 6, dpi = 300)





# ############# Functional analysis #############

## ########### Gene lists ###########

resLFC_df <- read.table("example_2_conditions/treatment_vs_control_result.tsv", 
                        sep = "\t", header = T)

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- log2(2)

# Background dataset (all genes tested for significance)
all_genes <- as.character(resLFC_df$gene)
str(all_genes)

# Significant DE genes
sig_genes <- as.character(
  resLFC_df %>%
  filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff) %>%
  pull(gene)
)

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
  gene = up_genes,        # gene list in SGD ORF IDs
  universe = all_genes,
  OrgDb = org.Sc.sgd.db,  # OrgDb for S. cerevisiae
  keyType = "ORF",        # Key type of gene identifiers
  ont = "BP",             # Ontology: BP (Biological Process)
  pAdjustMethod = "BH",   # Multiple testing correction
  pvalueCutoff  = 0.05,   # P-value cutoff
  qvalueCutoff  = 0.2     # Q-value cutoff
)

# Output results to a dataframe
ora_bp_up_df <- data.frame(ora_bp_up)
range(ora_bp_up_df$p.adjust, na.rm = TRUE)

# Save to file
write.table(ora_bp_up_df, 
            file = "example_2_conditions/treatment_vs_control_result_DE_up_ORA_BP.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Top n terms by p.adjust
# Display gene ratio (number of genes related to GO term / total number of significant genes)
dotplot(ora_bp_up, showCategory = 10, title="GO Biological Processes Enrichment in Up-regulated genes")
# Display gene count
png("example_2_conditions/treatment_vs_control_result_DE_up_ORA_BP_top10.png", width = 8, height = 6, units = "in", res = 300)
barplot(ora_bp_up, showCategory = 10, title="GO Biological Processes Enrichment in Up-regulated genes")
dev.off()

ora_bp_up_df %>%
  arrange(p.adjust) %>%
  head(n=10) %>%
  pull(Description)



### ######### Down regulated genes #########
ora_bp_down <- enrichGO(
  gene = down_genes,        # gene list in SGD ORF IDs
  universe = all_genes,
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
write.table(ora_bp_down_df, 
            file = "example_2_conditions/treatment_vs_control_result_DE_down_ORA_BP.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Top n terms by p.adjust
# Display gene ratio (number of genes related to GO term / total number of significant genes)
dotplot(ora_bp_down, showCategory = 10, title="GO Biological Processes Enrichment in Down-regulated genes")
# Display gene count
png("example_2_conditions/treatment_vs_control_result_DE_down_ORA_BP_top10.png", width = 8, height = 6, units = "in", res = 300)
barplot(ora_bp_down, showCategory = 10, title="GO Biological Processes Enrichment in Down-regulated genes")
dev.off()

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


