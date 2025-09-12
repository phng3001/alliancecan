#!/bin/bash
# P=NP

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: bash $0 <snpEff_container_path> <genome_name> <ref_fasta_path> <ref_gff_path> <mutant_list>"
    exit 1
fi

# Assign arguments to variables
snpEff_container="$1"
genome="$2"
ref_fasta="$3"
ref_gff="$4"
mutant_list="$5"

# Prerequisites
scripts=("aggregate_variant.py" \
"combine_variant.py" \
"calculate_allele_frequency.py" \
"add_gene_description.py" \
"compare_variant.py" \
"add_allele_frequency.py" \
"build_presence_absence_matrix.py")

# Check if the prerequisite scripts exist in the current directory
for script in "${scripts[@]}"
do
    if [ ! -f "$script" ]; then
    echo "Error: Required script $script not found in the current directory."
    exit 1
    fi
done

# Load modules
module purge
module load StdEnv/2023 apptainer/1.2.4 scipy-stack/2024a

# Declare variables
realpath=$(realpath .)
tmp_dir="tmp"

# Make temporary folder 
if [ -e "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir

# Clean old files if exist
shopt -s nullglob

files=(*_filtered_variants.noWT.snpEff.gatk.vcf \
*_filtered_variants.noWT.snpEff.freebayes.vcf \
*_filtered_variants.noWT.snpEff.bcftools.vcf)

if [ ${#files[@]} -gt 0 ]; then
    echo "Moving the following old files to $tmp_dir:"
    printf '%s\n' "${files[@]}"
    mv "${files[@]}" $tmp_dir
fi

# Handle gff file 
sed '/^$/d' $ref_gff > ${genome}.gff # remove empty lines (NCBI gff)
sed -i '/##FASTA/,$d' ${genome}.gff # remove fasta sequence (Prokka gff)



######### SnpEff annotation #########

# Build SnpEff database
echo "Start building SnpEff database"
echo -e "data.dir=/opt/snpEff_db/\n${genome}.genome:$genome" > snpEff.config
mkdir -p snpEff_db/$genome
cp $ref_fasta snpEff_db/$genome/sequences.fa
cp ${genome}.gff snpEff_db/$genome/genes.gff
apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
snpEff build -gff3 -v -noCheckCds -noCheckProtein $genome
echo "SnpEff database built"

# Annote mutant variants with SnpEff
## GATK variants
for X in $(cat $mutant_list)
do
    echo "Annote sample $X's GATK variants"
    apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
    snpEff ann -v -no-downstream -no-upstream $genome ${X}_filtered_variants.noWT.gatk.vcf \
    > ${X}_filtered_variants.noWT.snpEff.gatk.vcf
done

## FreeBayes variants
for X in $(cat $mutant_list)
do
    echo "Annote sample $X's FreeBayes variants"
    apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
    snpEff ann -v -no-downstream -no-upstream $genome ${X}_filtered_variants.noWT.freebayes.vcf \
    > ${X}_filtered_variants.noWT.snpEff.freebayes.vcf
done

## Bcftools variants
for X in $(cat $mutant_list)
do
    echo "Annote sample $X's Bcftools variants"
    apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
    snpEff ann -v -no-downstream -no-upstream $genome ${X}_filtered_variants.noWT.bcftools.vcf \
    > ${X}_filtered_variants.noWT.snpEff.bcftools.vcf
done



######### SnpEff postprocessing #########

# Aggregate variants from all samples
# Combine alternate alleles in multi-allelic sites
# Calculate allele frequencies
## GATK variants
echo "Start aggregating GATK variants"
python aggregate_variant.py *_filtered_variants.noWT.snpEff.gatk.vcf > all_variant.gatk.tmp0
sed -i 's/_filtered_variants.noWT.snpEff.gatk.vcf//g' all_variant.gatk.tmp0 
sed -i '/non_coding_transcript_variant/d' all_variant.gatk.tmp0 # pseudogenes
sed -i '/intragenic_variant/d' all_variant.gatk.tmp0 # redundant
python combine_variant.py -i all_variant.gatk.tmp0 -o all_variant.gatk.tmp1
python calculate_allele_frequency.py -i all_variant.gatk.tmp1 -o all_variant.gatk.tsv
rm all_variant.gatk.tmp0 all_variant.gatk.tmp1

## FreeBayes variants
echo "Start aggregating FreeBayes variants"
python aggregate_variant.py *_filtered_variants.noWT.snpEff.freebayes.vcf > all_variant.freebayes.tmp0
sed -i 's/_filtered_variants.noWT.snpEff.freebayes.vcf//g' all_variant.freebayes.tmp0
sed -i '/non_coding_transcript_variant/d' all_variant.freebayes.tmp0 # pseudogenes
sed -i '/intragenic_variant/d' all_variant.freebayes.tmp0 # redundant
python combine_variant.py -i all_variant.freebayes.tmp0 -o all_variant.freebayes.tmp1
python calculate_allele_frequency.py -i all_variant.freebayes.tmp1 -o all_variant.freebayes.tsv
rm all_variant.freebayes.tmp0 all_variant.freebayes.tmp1

## Bcftools variants
echo "Start aggregating Bcftools variants"
python aggregate_variant.py *_filtered_variants.noWT.snpEff.bcftools.vcf > all_variant.bcftools.tmp0
sed -i 's/_filtered_variants.noWT.snpEff.bcftools.vcf//g' all_variant.bcftools.tmp0
sed -i '/non_coding_transcript_variant/d' all_variant.bcftools.tmp0 # pseudogenes
sed -i '/intragenic_variant/d' all_variant.bcftools.tmp0 # redundant
python combine_variant.py -i all_variant.bcftools.tmp0 -o all_variant.bcftools.tmp1
python calculate_allele_frequency.py -i all_variant.bcftools.tmp1 -o all_variant.bcftools.tsv
rm all_variant.bcftools.tmp0 all_variant.bcftools.tmp1


# Add gene description
## GATK variants
echo "Start adding gene description for GATK variants"
cp all_variant.gatk.tsv all_variant.gatk.tmp
python add_gene_description.py \
  -i all_variant.gatk.tmp \
  -g snpEff_db/$genome/genes.gff \
  -o all_variant.gatk.tsv
rm all_variant.gatk.tmp

## FreeBayes variants
echo "Start adding gene description for FreeBayes variants"
cp all_variant.freebayes.tsv all_variant.freebayes.tmp
python add_gene_description.py \
  -i all_variant.freebayes.tmp \
  -g snpEff_db/$genome/genes.gff \
  -o all_variant.freebayes.tsv
rm all_variant.freebayes.tmp

## Bcftools variants
echo "Start adding gene description for Bcftools variants"
cp all_variant.bcftools.tsv all_variant.bcftools.tmp
python add_gene_description.py \
  -i all_variant.bcftools.tmp \
  -g snpEff_db/$genome/genes.gff \
  -o all_variant.bcftools.tsv
rm all_variant.bcftools.tmp


# Combine GATK, FreeBayes and Bcftools variant information
# Consider variants detected at the same position as consensus
echo "Start combining GATK, FreeBayes and Bcftools variant information"
key_columns="#CHROM POS GENOTYPE GENE_ID GENE_NAME DESCRIPTION FEATURE_ID VARIANT_TYPE NU_CHANGE AA_CHANGE CDS_POS/CDS_LENGTH AA_POS/AA_LENGTH SAMPLE"
python compare_variant.py \
  -k $key_columns \
  -n gatk -f all_variant.gatk.tsv \
  -n freebayes -f all_variant.freebayes.tsv \
  -n bcftools -f all_variant.bcftools.tsv \
  -o all_variant.tmp


# Add allele frequency information
# Order of priority: GATK > FreeBayes > Bcftools
key_columns="#CHROM POS NU_CHANGE SAMPLE"
python add_allele_frequency.py \
  -i all_variant.tmp \
  -r all_variant.gatk.tsv all_variant.freebayes.tsv all_variant.bcftools.tsv \
  -k $key_columns \
  -c ALLELE_FREQ \
  -p 3 \
  -o all_variant.tsv \
  -v
rm all_variant.tmp


# Build presence/absence matrix across samples
echo "Start building mutated gene presence/absence matrix"
key_columns="#CHROM GENE_ID GENOTYPE GENE_NAME DESCRIPTION VARIANT_TYPE CALLED_BY"
python build_presence_absence_matrix.py \
  -i all_variant.tsv \
  -k $key_columns \
  -s SAMPLE \
  -o all_variant_matrix.tsv


