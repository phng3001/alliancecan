#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: bash $0 <snpEff_container_path> <genome_name> <ref_fasta_path> <ref_gff_path> <sample_list>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.4.5

# Assign arguments to variables
snpEff_container="$1"
genome="$2"
ref_fasta="$3"
ref_gff="$4"
sample_list="$5"

# Declare variables
realpath=$(realpath .)
extension=".filt.vcf"

# Handle gff file 
sed '/^$/d' $ref_gff > ${genome}.gff # remove empty lines (NCBI gff)
sed -i '/##FASTA/,$d' ${genome}.gff # remove fasta sequence (Prokka gff)

# Build SnpEff database
echo "Start building SnpEff database"
echo -e "data.dir=/opt/snpEff_db/\n${genome}.genome:$genome" > snpEff.config
mkdir -p snpEff_db/$genome
cp $ref_fasta snpEff_db/$genome/sequences.fa
cp ${genome}.gff snpEff_db/$genome/genes.gff
apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
snpEff build -gff3 -v -noCheckCds -noCheckProtein $genome
echo "SnpEff database built"

# Annote variants with SnpEff
for X in $(cat $sample_list)
do
    echo "Annote sample $X variants"
    apptainer exec -B $realpath/snpEff_db:/opt/snpEff_db $snpEff_container \
    snpEff ann -v -no-downstream -no-upstream $genome ${X}${extension} \
    > ${X}.filt.snpeff.vcf
    mv snpEff_genes.txt ${X}_snpEff_genes.txt
    mv snpEff_summary.html ${X}_snpEff_summary.html
done


