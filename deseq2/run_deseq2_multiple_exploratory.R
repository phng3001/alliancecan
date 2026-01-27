# P=NP

cat("Hello from P=NP!\n")
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_name <- basename(script_path)
cat("Running script:", script_name, "\n")

# ----- Parse trailing args -----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript ", script_name, " <count_data> <col_data> <design> [transcriptId_col] [sampleId_col] [condition_col] [transform_function]")
}

# Required arguments
count_data <- args[1]
col_data <- args[2]
design_arg <- args[3]

# Optional arguments
transcriptId_col <- ifelse(length(args) >= 4, args[4], "transcript_id")
sampleId_col <- ifelse(length(args) >= 5, args[5], "id")
condition_col <- ifelse(length(args) >= 6, args[6], "condition")
transform_function <- ifelse(length(args) >= 7, args[7], "vst")

# Add ~ if not present in design
if (!grepl("~", design_arg)) {
  design_arg <- paste0("~", design_arg)
}
design_formula <- as.formula(design_arg)

cat("Arguments:\n")
cat("  Count data:", count_data, "\n")
cat("  Col data:", col_data, "\n")
cat("  Design:", design_arg, "\n")
cat("  Transcript column:", transcriptId_col, "\n")
cat("  Sample ID column:", sampleId_col, "\n")
cat("  Condition column:", condition_col, "\n")
cat("  Transformation:", transform_function, "\n\n")



# ----- Load packages -----
packages <- c("DESeq2", "pheatmap", "RColorBrewer")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package failed to load:", pkg))
  }
}



# ----- Read input data -----
# countData <- read.csv(count_data, row.names = transcriptId_col)
# colData <- read.csv(col_data, row.names = sampleId_col)

countData <- read.csv(count_data)
colData <- read.csv(col_data)

# Check that transcript column exists
if (!(transcriptId_col %in% colnames(countData))) {
  stop(paste("Transcript ID column", transcriptId_col, "not found in count data file"))
}

# Check that sample ID column exists
if (!(sampleId_col %in% colnames(colData))) {
  stop(paste("Sample ID column", sampleId_col, "not found in col data file"))
}

# Set row names
rownames(countData) <- countData[[transcriptId_col]]
countData[[transcriptId_col]] <- NULL

rownames(colData) <- colData[[sampleId_col]]
colData[[sampleId_col]] <- NULL

# Check sample ID match
if (!all(colnames(countData) == rownames(colData))) {
  stop("Column names in count data file do not match row names in col data file")
}



# ---- Build DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = design_formula)
cat("Raw DESeqDataSet:\n")
print(dds)
cat("\n")

# Pre-filter: remove low count genes
# keep <- rowSums(counts(dds)) >= 10
# dds <- dds[keep,]
# cat("Filtered DESeqDataSet (min 10 reads):\n")
# print(dds)
# cat("\n")

cond <- colData(dds)[[condition_col]]
keep <- rep(FALSE, nrow(dds))
for (c in levels(cond)) {
  samples_in_c <- which(cond == c)
  keep <- keep | rowSums(counts(dds)[, samples_in_c, drop=FALSE]) >= 10
}
dds <- dds[keep,]
cat("Filtered DESeqDataSet (min 10 reads in at least one condition):\n")
print(dds)
cat("\n")



# ---- Run DESeq ----
dds <- DESeq(dds)
res <- results(dds)



# ----- Transform counts -----
n_genes <- nrow(dds)
cat("Number of genes after filtering:", n_genes, "\n")

# Count how many genes have mean normalized count > 5
expressed_genes <- sum(rowMeans(counts(dds, normalized = TRUE)) > 5)

if (transform_function == "vst") {
  if (expressed_genes < 1000) {
    cat("Using varianceStabilizingTransformation (fewer than 1000 expressed genes with mean normalized count > 5)\n")
    trans <- varianceStabilizingTransformation(dds, blind = TRUE)
  } else {
    cat("Using vst()\n")
    trans <- vst(dds, blind = TRUE)
  }
} else if (transform_function == "rlog") {
  cat("Using rlog()\n")
  trans <- rlog(dds, blind = TRUE)
} else {
  stop("Invalid transformation method (must be 'vst' or 'rlog')")
}



# ----- Visualization -----
# PCA plot
plot_name <- "PCA_plot.pdf"
pdf(plot_name)
plotPCA(trans, intgroup = all.vars(design_formula)[1])
dev.off()
cat("PCA plot saved as", plot_name, "\n")

# Sample distance heatmap
plot_name <- "heatmap_sample_distance.pdf"
pdf(plot_name)
sampleDists <- dist(t(assay(trans)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(trans[[all.vars(design_formula)[1]]])
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
cat("Sample-to-sample distance heatmap saved as", plot_name, "\n")