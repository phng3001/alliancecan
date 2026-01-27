#!/bin/bash
# P=NP

# Prerequisites
scripts=("run_deseq2.R")

# Check if the prerequisite scripts exist in the current directory
for script in "${scripts[@]}"
do
    if [ ! -f "$script" ]; then
    echo "Error: Required script $script not found in the current directory."
    exit 1
    fi
done

# Check if the number of arguments is correct
# Number of mandatory and optional arguments
MANDATORY_ARGS=7
TOTAL_ARGS=13

if [ $# -lt $MANDATORY_ARGS ]; then
    echo "Error: You must provide at least $MANDATORY_ARGS arguments."
    echo "Usage: bash $0 <count_data> <col_data> <design> <ref_variable> <ref_condition> <transcript_annotation> <output_prefix> [transcriptId_col] [sampleId_col] [condition_col] [top_n] [transform_function] [shrinkage_coef]"
    exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
    echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."        
    echo "Usage: bash $0 <count_data> <col_data> <design> <ref_variable> <ref_condition> <transcript_annotation> <output_prefix> [transcriptId_col] [sampleId_col] [condition_col] [top_n] [transform_function] [shrinkage_coef]"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 r/4.4.0 r-bundle-bioconductor/3.20

# Asign arguments to variables
count_data="$1"
col_data="$2"
design="$3"
ref_variable="$4"
ref_condition="$5"
transcript_annotation="$6"
output_prefix="$7"
transcriptId_col="${8:-transcript_id}"
sampleId_col="${9:-id}"
condition_col="${10:-condition}"
top_n="${11:-30}"
transform_function="${12:-vst}"
shrinkage_coef="${13}"

# Run R script
Rscript \
run_deseq2.R \
$count_data \
$col_data \
$design \
$ref_variable \
$ref_condition \
$transcript_annotation \
$output_prefix \
$transcriptId_col \
$sampleId_col \
$condition_col \
$top_n \
$transform_function \
$shrinkage_coef