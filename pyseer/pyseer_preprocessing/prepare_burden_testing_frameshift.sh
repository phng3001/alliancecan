#!/bin/bash
# P=NP
# VCF files annotated with SnpEff are required for this script to work properly

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



######### Variant filtering #########

# Filter VCF for frameshift variants
bcftools view \
-i 'ANN[*] ~ "frameshift_variant"' "$merged_vcf" \
| bgzip -c > ${output_prefix}.tmp0.vcf.gz

tabix -p vcf ${output_prefix}.tmp0.vcf.gz # indexing

# Fix the genotype field from ./. (absence) to 0/0 (reference)
bcftools +setGT \
${output_prefix}.tmp0.vcf.gz \
-Oz -o ${output_prefix}.tmp1.vcf.gz \
-- -t q -n 0 \
-i 'GT="./."'

tabix -p vcf ${output_prefix}.tmp1.vcf.gz # indexing

# Only keep GT among the FORMAT fields
bcftools annotate \
-x FORMAT/DP,FORMAT/RO,FORMAT/QR,FORMAT/AO,FORMAT/QA,FORMAT/GL \
${output_prefix}.tmp1.vcf.gz \
-Oz -o ${output_prefix}.tmp2.vcf.gz

tabix -p vcf ${output_prefix}.tmp2.vcf.gz # indexing

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
${output_prefix}.tmp0.vcf.gz.tbi \
${output_prefix}.tmp1.vcf.gz \
${output_prefix}.tmp1.vcf.gz.tbi \
${output_prefix}.tmp2.vcf.gz \
${output_prefix}.tmp2.vcf.gz.tbi

if [ -s "${output_prefix}.vcf.gz" ]; then
    echo "Frameshift variant filtering completed. Output files:"
    echo "${output_prefix}.vcf.gz"
    echo "${output_prefix}.vcf.gz.tbi"
else
    echo "Problem filtering frameshift variants, output file ${output_prefix}.vcf.gz was not generated."
fi



######### Variant region extraction #########

# Get variant gene list
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' \
${output_prefix}.vcf.gz | \
awk -F'\t' '
{
    n = split($3, ann, ",")

    for (i=1; i<=n; i++) {
        split(ann[i], field, "|")

        annotation = field[2]
        locus_tag = field[5]

        if (annotation ~ /frameshift_variant/) {
            print locus_tag "\t" $1 ":" $2 "-" $2
        }
    }
}' > ${output_prefix}_variant_positions.txt

cut -f1 ${output_prefix}_variant_positions.txt \
| sort | uniq > ${output_prefix}_variant_genes.txt

# Get gene regions from GFF
grep -F -f ${output_prefix}_variant_genes.txt "$reference_gff" \
> ${output_prefix}_variant_genes.gff

awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene" {
    match($9,/ID=([^;]+)/,a)
    if(a[1]!="")
        print a[1], $1":"$4"-"$5
}
' ${output_prefix}_variant_genes.gff > ${output_prefix}_variant_gene_regions.txt

# Pad gene regions by 5 bp on each side as variants may occur at the edges of genes
awk -F'\t' '
BEGIN{OFS="\t"}
$3=="gene" {
    match($9,/ID=([^;]+)/,a)
    if(a[1]!="") {
        start = ($4 > 5 ? $4 - 5 : 1)
        end = $5 + 5
        print a[1], $1 ":" start "-" end
    }
}
' ${output_prefix}_variant_genes.gff > ${output_prefix}_variant_gene_regions_padded_5bp.txt

# Remove temporary files
rm \
${output_prefix}_variant_positions.txt \
${output_prefix}_variant_genes.txt \
${output_prefix}_variant_genes.gff 

if [[ -s "${output_prefix}_variant_gene_regions.txt" && -s "${output_prefix}_variant_gene_regions_padded_5bp.txt" ]]; then
    echo "Variant region extraction completed. Output files:"
    echo "${output_prefix}_variant_gene_regions.txt"
    echo "${output_prefix}_variant_gene_regions_padded_5bp.txt"
else
    echo "Problem extracting variant regions, output files were not generated."
fi
