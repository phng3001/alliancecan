# P=NP

cat("Hello from P=NP!\n")
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_name <- basename(script_path)
cat("Running script:", script_name, "\n")

# ----- Parse trailing args -----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("Usage: Rscript ", script_name, " <count_data> <col_data> <design> <ref_variable> <ref_condition> <transcript_annotation> <output_prefix> [transcriptId_col] [sampleId_col] [top_n] [transform_function] [shrinkage_coef]")
}

# Required arguments
count_data <- args[1]
col_data <- args[2]
design_arg <- args[3]
ref_variable <- args[4]
ref_condition <- args[5]
transcript_annotation <- args[6]
output_prefix <- args[7]

# Optional arguments
transcriptId_col <- ifelse(length(args) >= 8, args[8], "transcript_id")
sampleId_col <- ifelse(length(args) >= 9, args[9], "id")
condition_col <- ifelse(length(args) >= 10, args[10], "condition")
top_n <- ifelse(length(args) >= 11, args[11], 30)
transform_function <- ifelse(length(args) >= 12, args[12], "vst")
shrinkage_coef <- ifelse(length(args) >= 13, args[13], NA)

# Add ~ if not present in design
if (!grepl("~", design_arg)) {
  design_arg <- paste0("~", design_arg)
}
design_formula <- as.formula(design_arg)

cat("Arguments:\n")
cat("  Count data:", count_data, "\n")
cat("  Col data:", col_data, "\n")
cat("  Design:", design_arg, "\n")
cat("  Reference variable:", ref_variable, "\n")
cat("  Reference condition:", ref_condition, "\n")
cat("  Transcript column:", transcriptId_col, "\n")
cat("  Sample ID column:", sampleId_col, "\n")
cat("  Condition column:", condition_col, "\n")
cat("  Transformation:", transform_function, "\n\n")
cat("  Number of top genes to be plotted:", top_n, "\n\n")


# ----- Load packages -----
packages <- c("DESeq2", "pheatmap", "RColorBrewer", "ReportingTools")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package failed to load:", pkg))
  }
}

# ----- Read input data -----
countData <- read.csv(count_data)
colData <- read.csv(col_data)
annotData <- read.csv(transcript_annotation)

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

# Check transcript annotation data
required_cols <- c(transcriptId_col, "chrom", "size", "description")
missing_cols <- setdiff(required_cols, colnames(annotData))
if (length(missing_cols) > 0) {
  stop("The following required columns are missing in the annotation file: ", paste(missing_cols, collapse = ", "))
}

# ----- Build DESeq2 dataset -----
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

# Set reference level
design_vars <- all.vars(design_formula)
if (!(ref_variable %in% design_vars)) {
  stop(paste("Reference variable", ref_variable, "not found in design formula"))
}
if (!(ref_condition %in% colData[[ref_variable]])) {
  stop(paste("Reference condition", ref_condition, "not found in col data column", ref_variable))
}
dds[[ref_variable]] <- relevel(dds[[ref_variable]], ref = ref_condition)
cat("Set reference level for", ref_variable, "to", ref_condition, "\n")



# ----- Run DESeq -----
dds <- DESeq(dds)
res <- results(dds)
res[[transcriptId_col]] <- rownames(res)

# Add annotation
annotData <- annotData[, c("transcript_id", "chrom", "size", "description")]
annotData_align <- annotData[match(res[[transcriptId_col]], annotData[[transcriptId_col]]), ]
res$chrom <- annotData_align$chrom
res$size <- annotData_align$size
res$description <- annotData_align$description
res[[transcriptId_col]] <- NULL # now redundant

# Coefficient selection for shrinkage
coef_list <- resultsNames(dds)
cat("Available coefficients:\n")
print(coef_list)

if (!is.na(shrinkage_coef)) {
  if (!(shrinkage_coef %in% coef_list)) {
    stop(paste("Provided shrinkage_coef", shrinkage_coef, "is not an available coefficient"))
  }
  shrink_name <- shrinkage_coef
  cat("Using user-supplied coefficient:", shrink_name, "\n")
} else {
  shrink_name <- coef_list[grepl(paste0(ref_variable, "_"), coef_list)]
  
  if (length(shrink_name) == 0) {
    cat("\nCould not auto-detect shrinkage coefficient.\n")
    cat("Available coefficients:\n")
    for (i in seq_along(coef_list)) {
      cat(sprintf("  [%d] %s\n", i, coef_list[i]))
  }

  if (!interactive()) {
    stop("Non-interactive session: Please rerun with the correct coefficient passed as a 10th argument.")
  }

  selection <- as.integer(readline(prompt = "Please enter the number of the coefficient to shrink: "))
  if (is.na(selection) || selection < 1 || selection > length(coef_list)) {
    stop("Invalid selection. Aborting.")
  }
  shrink_name <- coef_list[selection]
  }
  cat("Automatically selected shrinkage coefficient:", shrink_name, "\n")
}  

# LFC shrinkage
resLFC <- lfcShrink(dds, coef = shrink_name, type = "apeglm")



# ----- Summarize results -----
# res <- results(dds)
# res0.05 <- results(dds, alpha = 0.05)
# summary(res0.05)
summary(res)

# Sort by padj 
resOrdered <- res[order(res$padj), ]
write.csv(as.data.frame(resOrdered), paste0(output_prefix, "_result.csv"))

# Filter for DE genes
# FDR < 0.05 & |LFC| > 1
# resDE <- resOrdered[resOrdered$padj < 0.05 & abs(resOrdered$log2FoldChange) > 1, ]
resDE <- resOrdered[!is.na(resOrdered$padj) &
                    resOrdered$padj < 0.05 &
                    abs(resOrdered$log2FoldChange) > 1, ]

cat("DE genes (FDR < 0.05 & |LFC| > 1):\n")
print(resDE)
cat("\n")
write.csv(as.data.frame(resDE), paste0(output_prefix, "_result_DE.csv"))



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
# MA plot
plot_name <- paste0(output_prefix, "_MA_plot.pdf")
pdf(plot_name)
plotMA(resLFC, ylim = c(-5, 5))
dev.off()
cat("MA plot saved as", plot_name, "\n")

# PCA plot
plot_name <- paste0(output_prefix, "_PCA_plot.pdf")
pdf(plot_name)
plotPCA(trans, intgroup = ref_variable)
dev.off()
cat("PCA plot saved as", plot_name, "\n")

# Sample distance heatmap
plot_name <- paste0(output_prefix, "_heatmap_sample_distance.pdf")
pdf(plot_name)
sampleDists <- dist(t(assay(trans)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colData[[ref_variable]])
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
dev.off()
cat("Sample-to-sample distance heatmap saved as", plot_name, "\n")

# Top expressed heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:top_n]
plot_name <- paste0(output_prefix, "_heatmap_top", top_n, "_expressed.pdf")
pdf(plot_name)
pheatmap(assay(trans)[select, ],
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         annotation_col = colData)
dev.off()
cat("Transformed count matrix heatmap of top", top_n, "expressed genes saved as", plot_name, "\n")

# Top DE heatmap
select <- order(res$padj)[1:top_n]
plot_name <- paste0(output_prefix, "_heatmap_top", top_n,"_DE.pdf")
pdf(plot_name)
pheatmap(assay(trans)[select, ],
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE,
         annotation_col = colData)
dev.off()
cat("Transformed count matrix heatmap of top", top_n ,"DE genes saved as", plot_name, "\n")



# ----- ReportingTools -----
# Define the output directory name
dir_name <- paste0(output_prefix, "_reportingtools")
# Check if the directory exists
if (dir.exists(dir_name)) {
  # Remove the existing directory and its contents
  unlink(dir_name, recursive = TRUE)
  cat("Directory", dir_name, "has been reinitialized", "\n")
}
# Create the directory
dir.create(dir_name)
# Get the directory path
path <- file.path(getwd(), dir_name)

cat("Start reporting DE genes:\n")
# DE genes
## Annotation report
report_annot <- HTMLReport(shortName="DE_genes",
                           title=paste0("DE Genes: ", output_prefix),
                           basePath=path,
                           reportDirectory="FDR_0.05_LFC_1_annot")
if (nrow(resDE) > 0) {
  resDE_df <- as.data.frame(resDE)
  resDE_df[[transcriptId_col]] <- rownames(resDE)
  rownames(resDE_df) <- NULL
  resDE_df <- resDE_df[, c(transcriptId_col, setdiff(names(resDE_df), transcriptId_col))]
  publish(resDE_df, report_annot)
} else {
  cat("No DE genes passed the threshold; skipping annotation report.\n")
}
finish(report_annot)

## Plot report
report_plot <- HTMLReport(shortName="DE_genes",
                          title=paste0("DE Genes: ", output_prefix),
                          basePath=path,
                          reportDirectory="FDR_0.05_LFC_1_plot")
if (nrow(resDE) > 0) {
  tryCatch({
    publish(resDE, report_plot, dds, make.plots=TRUE, factor=dds[[ref_variable]])
    finish(report_plot)
  }, error = function(e) {
    cat("Skipping plot report for DE genes: No features meet selection criteria for plotting.\n")
  })
} else {
  cat("No DE genes passed the threshold; skipping plot report.\n")
}

cat("Start reporting top", top_n, "genes:\n")
# Top genes
resOrdered_df <- as.data.frame(resOrdered)
resOrdered_df[[transcriptId_col]] <- rownames(resOrdered)
rownames(resOrdered_df) <- NULL
resOrdered_df <- resOrdered_df[, c(transcriptId_col, setdiff(names(resOrdered_df), transcriptId_col))]

actual_top_n <- min(as.numeric(top_n), nrow(resOrdered_df))
if (actual_top_n < as.numeric(top_n)) {
  cat(paste0("Requested top_n (", top_n, ") exceeds number of available results (", nrow(resOrdered_df), "). Showing top ", actual_top_n, " instead.\n"))
}

## Annotation report
if (nrow(resOrdered) > 0) {
  report_annot <- HTMLReport(shortName=paste0("top_", actual_top_n, "_genes"),
                             title=paste0("Top ", actual_top_n, " genes:", output_prefix),
                             basePath=path,
                             reportDirectory=paste0("top_", actual_top_n, "_annot"))
  top_res <- head(resOrdered_df, n = actual_top_n)
  publish(top_res, report_annot)
  finish(report_annot)
} else {
  cat("Skipping annotaion report: No results available.\n")
}

## Plot report
if (nrow(resOrdered) > 0) {
  report_plot <- HTMLReport(shortName=paste0("top_", actual_top_n, "_genes"),
                            title=paste0("Top ", actual_top_n, " genes:", output_prefix),
                            basePath=path,
                            reportDirectory=paste0("top_", actual_top_n, "_plot"))
  
  tryCatch({
    publish(res, report_plot, dds, n=top_n, make.plots=TRUE, factor=dds[[ref_variable]])
    finish(report_plot)
  }, error = function(e) {
    cat("Skipping plot report: No features meet selection criteria for plotting.\n")
  })
} else {
  cat("Skipping plot report: No results available.\n")
}
