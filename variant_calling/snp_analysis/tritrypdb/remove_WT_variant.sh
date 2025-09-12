#!/bin/bash
# P=NP

# Check if the remove_WT_variant.py script exists in the current directory
if [ ! -f "remove_WT_variant.py" ]; then
    echo "Error: Required script remove_WT_variant.py not found in the current directory."
    exit 1
fi

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: bash $0 <calling_folder> <wt_sample_name> <mutant_list>"
    exit 1
fi

# Assign arguments to variables
var_calling_folder="$1"
wt="$2"
mutant_list="$3"

# Remove WT variants from mutant variant lists
# gatk variants
while read -r X; do
    echo "Remove WT gatk variants & N-stretch-containing variants from sample $X"
    python remove_WT_variant.py \
    -w "$var_calling_folder/$wt/${wt}_filtered_variants.gatk.vcf" \
    -m "$var_calling_folder/$X/${X}_filtered_variants.gatk.vcf" \
    -o "${X}_filtered_variants.noWT.gatk.vcf"
done < "$mutant_list"

# freebayes variants
while read -r X; do
    echo "Remove WT freebayes variants & N-stretch-containing variants from sample $X"
    python remove_WT_variant.py \
    -w "$var_calling_folder/$wt/${wt}_filtered_variants.freebayes.vcf" \
    -m "$var_calling_folder/$X/${X}_filtered_variants.freebayes.vcf" \
    -o "${X}_filtered_variants.noWT.freebayes.vcf"
done < "$mutant_list"

# bcftools variants
while read -r X; do
    echo "Remove WT bcftools variants & N-stretch-containing variants from sample $X"
    python remove_WT_variant.py \
    -w "$var_calling_folder/$wt/${wt}_filtered_variants.bcftools.vcf" \
    -m "$var_calling_folder/$X/${X}_filtered_variants.bcftools.vcf" \
    -o "${X}_filtered_variants.noWT.bcftools.vcf"
done < "$mutant_list"
