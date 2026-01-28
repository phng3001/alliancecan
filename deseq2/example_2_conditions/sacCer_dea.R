############# Initial setup #############

# Set working directory
setwd("/project/def-mouellet/Scripts_MOU/PNP/alliancecan/deseq2/")

# Activate environment
renv::activate()

# Load required packages
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)





############# Load data #############

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





############# Differential expression analysis #############

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = meta,
  design = ~condition
)

# Original count matrix
# counts(dds)
# str(counts(dds))

########### Pre-filter (optional) ###########
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



########### Run DESeq2 ###########

# Set reference level
dds$condition <- relevel(dds$condition, ref = "C")

# Run DE analysis
dds <- DESeq(dds)



########### Check the model fit ###########

# Dispersion (α): Var = μ + α*μ^2
# Dispersion ~ 1/μ, ~ Var
png("example_2_conditions/dispersion_estimates.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds)
dev.off()



########### Explore the results ###########

# Available coefficients
resultsNames(dds)

######### Unshrunken LFC #########

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



######### Shrunken LFC #########
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


