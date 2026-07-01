#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: bash $0 <ref_gff_path> <gtf_merge_list> <gap_len> <output_gtf> <gffcompare_outprefix>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 stringtie/3.0.1 gffcompare/0.12.6

# Asign arguments to variables
ref_gff_path="$1"
gtf_merge_list="$2"
gap_len="$3"
output_gtf="$4"
gffcompare_outprefix="$5"

# StringTie merge
## merge/assemble transcripts into a non-redundant set of transcripts
stringtie --merge \
-G $ref_gff_path \
-o $output_gtf \
-g $gap_len \
$gtf_merge_list

# gffcompare
## evaluate transcript assemblies
gffcompare \
-r $ref_gff_path \
-o $gffcompare_outprefix \
$output_gtf
