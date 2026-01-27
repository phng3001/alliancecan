#!/bin/bash
# P=NP

# Prerequisites
scripts=("run_deseq2_multiple_exploratory.R")

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
MANDATORY_ARGS=3
TOTAL_ARGS=7

if [ $# -lt $MANDATORY_ARGS ]; then
    echo "Error: You must provide at least $MANDATORY_ARGS arguments."
    echo "Usage on terminal: bash $0 <count_data> <col_data> <design> [transcriptId_col] [sampleId_col] [condition_col] [transform_function]"
    exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
    echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."        
    echo "Usage on terminal: bash $0 <count_data> <col_data> <design> [transcriptId_col] [sampleId_col] [condition_col] [transform_function]"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 r/4.4.0 r-bundle-bioconductor/3.20

# Asign arguments to variables
count_data="$1"
col_data="$2"
design="$3"
transcriptId_col="${4:-transcript_id}"
sampleId_col="${5:-id}"
condition_col="${6:-condition}"
transform_function="${7:-vst}"

# Run R script
Rscript \
run_deseq2_multiple_exploratory.R \
$count_data \
$col_data \
$design \
$transcriptId_col \
$sampleId_col \
$condition_col \
$transform_function