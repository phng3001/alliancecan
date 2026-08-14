#!/bin/bash
# P=NP
# Prerequisites
## Merged (multi-sample) VCF file annotated with SnpEff
## NCBI GFF3 file

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <merged_vcf> <reference_fasta> <reference_gff> <output_prefix>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 bcftools/1.22

# Assign arguments to variables
merged_vcf="$1"
reference_fasta="$2"
reference_gff="$3"
output_prefix="$4"

# Prerequisites scripts
scripts=("get_promoter_region_from_gff.py" \
)
# Check if the prerequisite scripts exist in the working directory
for script in "${scripts[@]}"
do
    if [ ! -f "$script" ]; then
    echo "Error: Required script $script not found in the current directory"
    exit 1
    fi
done

echo "Preparing promoter region variant burden testing files for pyseer..."



######### Variant filtering #########

# Get promoter regions from GFF file
python3 get_promoter_region_from_gff.py \
--gff $reference_gff \
--output ${output_prefix}_regions.bed

# Filter VCF for variants in promoter regions
bcftools view \
-R ${output_prefix}_regions.bed "$merged_vcf" \
-Oz -o ${output_prefix}.tmp0.vcf.gz

# Fix the genotype field from ./. (absence) to 0/0 (reference)
bcftools +setGT \
${output_prefix}.tmp0.vcf.gz \
-Oz -o ${output_prefix}.tmp1.vcf.gz \
-- -t q -n 0 \
-i 'GT="./."'

# Only keep GT among the FORMAT fields
bcftools annotate \
-x FORMAT/DP,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
${output_prefix}.tmp1.vcf.gz \
-Oz -o ${output_prefix}.tmp2.vcf.gz

# Split multiallelic variants into multiple biallelic variants
# E.g. A/C,G -> A/C and A/G
bcftools norm \
-m -both \
-f "$reference_fasta" \
${output_prefix}.tmp2.vcf.gz \
-Oz -o ${output_prefix}.vcf.gz

tabix -p vcf ${output_prefix}.vcf.gz # indexing

# Remove temporary files
rm \
${output_prefix}.tmp0.vcf.gz \
${output_prefix}.tmp1.vcf.gz \
${output_prefix}.tmp2.vcf.gz

if [ -s "${output_prefix}.vcf.gz" ]; then
    echo "Promoter region variant filtering completed. Output files:"
    echo "${output_prefix}.vcf.gz"
    echo "${output_prefix}.vcf.gz.tbi"
else
    echo "Problem filtering promoter region variants, output file ${output_prefix}.vcf.gz was not generated."
fi



######### Variant region extraction #########

awk -F'\t' '{print $4, $1 ":" $2 "-" $3}' \
${output_prefix}_regions.bed \
> ${output_prefix}_regions.txt

if [ -s "${output_prefix}_regions.txt" ]; then
    echo "Promoter variant region extraction completed. Output file:"
    echo "${output_prefix}_regions.txt"
else
    echo "Problem extracting promoter region variant regions, output file ${output_prefix}_regions.txt was not generated."
fi
