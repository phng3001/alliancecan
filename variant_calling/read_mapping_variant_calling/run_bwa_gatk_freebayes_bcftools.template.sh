#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___bwa_gatk_freebayes_bcftools

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <ref_fasta_path> <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <ref_fasta_path> <fastq_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 samtools/1.20 bwa/0.7.18 gatk/4.4.0.0 r/4.4.0 bcftools/1.19 freebayes/1.3.7

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}
export spark_memory=4g
export spark_memory_overhead=2g # by experience

# Prerequisites
scripts=("combine_vcf_columns.py")

# Check if the prerequisite scripts exist in the current directory
for script in "${scripts[@]}"
do
    if [ ! -f "$script" ]; then
    echo "Error: $script not found in the current directory"
    exit 1
    fi
done

# Asign arguments to variables
ref_fasta_path="$1"
fastq_dir_path="$2"

# Declare variables
output_dir=__SAMPLE__
R1_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
R2_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)
tmp_dir=$HOME/scratch/${output_dir}_bwa_gatk_freebayes_bcftools

# Make temporary directory
## Check if the temporary directory exists and reset it
if [ -d "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir

# Make output directory
## Check if the output directory exists, if yes move it to the temporary directory
if [ -d "$output_dir" ]; then
    mv "$output_dir" "$tmp_dir"
    echo "Pre-existed folder $output_dir has been moved to $tmp_dir"
fi
mkdir $output_dir
mkdir $output_dir/reference
mkdir $output_dir/tmp

# Check if the fastq file paths are correct
fastq_paths=("$R1_path" "$R2_path")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done



######### Read mapping #########

# Index reference
## Get reference fasta file
cp $ref_fasta_path $output_dir/reference
## Get ref fasta file extension
ref_fasta_extension="${ref_fasta_path##*.}"
## Get ref fasta file base name without extension
ref_fasta_basename="${ref_fasta_path##*/}"
ref_fasta_basename="${ref_fasta_basename%.*}"
## BWA indexing
bwa index $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension}
## Samtools indexing
samtools faidx $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension}
## GATK indexing
gatk CreateSequenceDictionary R=$output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} O=$output_dir/reference/${ref_fasta_basename}.dict

# BWA alignment
bwa mem \
    -t $OMP_NUM_THREADS \
    -R "@RG\tID:__SAMPLE__\tPL:ILLUMINA\tSM:__SAMPLE__" \
    $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    $R1_path \
    $R2_path \
> $output_dir/__SAMPLE___paired.sam

# Check SAM file
samtools quickcheck $output_dir/__SAMPLE___paired.sam # no output means ok
gatk ValidateSamFile \
    -I $output_dir/__SAMPLE___paired.sam \
    -MODE VERBOSE

# Mark duplicates and sort
gatk MarkDuplicatesSpark \
    -I $output_dir/__SAMPLE___paired.sam \
    -O $output_dir/__SAMPLE___sorted_dedup_reads.bam \
    --spark-master local[$OMP_NUM_THREADS] \
    --conf spark.executor.memory=$spark_memory \
    --conf spark.executor.memoryOverhead=$spark_memory_overhead \
    --conf spark.executor.cores=$OMP_NUM_THREADS \
    --tmp-dir $TMPDIR

mv $output_dir/__SAMPLE___paired.sam $output_dir/tmp

# Collect alignment & insert size metrics
## Alignment metrics
gatk CollectAlignmentSummaryMetrics \
    R=$output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    I=$output_dir/__SAMPLE___sorted_dedup_reads.bam \
    O=$output_dir/__SAMPLE___alignment_metrics.txt
## Insert size metrics
gatk CollectInsertSizeMetrics \
    INPUT=$output_dir/__SAMPLE___sorted_dedup_reads.bam \
    OUTPUT=$output_dir/__SAMPLE___insert_size_metrics.txt \
    HISTOGRAM_FILE=$output_dir/__SAMPLE___insert_size_histogram.pdf



######### Variant calling #########

# --------------------
# GATK HaplotypeCaller
# --------------------

# Call variants
gatk HaplotypeCaller \
    -R $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    -I $output_dir/__SAMPLE___sorted_dedup_reads.bam \
    -O $output_dir/__SAMPLE___raw_variants.gatk.vcf

# Extract SNPs & INDELs
## SNPs
gatk SelectVariants \
    -R $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    -V $output_dir/__SAMPLE___raw_variants.gatk.vcf\
    --select-type SNP \
    -O $output_dir/__SAMPLE___raw_snps.gatk.vcf
## INDELs
gatk SelectVariants \
    -R $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    -V $output_dir/__SAMPLE___raw_variants.gatk.vcf\
    --select-type INDEL \
    -O $output_dir/__SAMPLE___raw_indels.gatk.vcf

# Filter variants
# https://pmc.ncbi.nlm.nih.gov/articles/PMC9163752/
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
## SNPs
gatk VariantFiltration \
    -R $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    -V $output_dir/__SAMPLE___raw_snps.gatk.vcf \
    -O $output_dir/__SAMPLE___filter_snps.gatk.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 60.0" \
    -filter-name "SOR_filter" -filter "SOR > 3.0" \
    -filter-name "MQ_filter" -filter "MQ < 40.0" \
    -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
    -genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 20"
## INDELs
gatk VariantFiltration \
    -R $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    -V $output_dir/__SAMPLE___raw_indels.gatk.vcf \
    -O $output_dir/__SAMPLE___filter_indels.gatk.vcf \
    -filter-name "QD_filter" -filter "QD < 2.0" \
    -filter-name "FS_filter" -filter "FS > 200.0" \
    -filter-name "SOR_filter" -filter "SOR > 10.0" \
    -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -20.0" \
    -filter-name "InbreedingCoeff_filter" -filter "InbreedingCoeff < -0.8" \
    -genotype-filter-name "DP_filter" -genotype-filter-expression "DP < 10" \
    -genotype-filter-name "GQ_filter" -genotype-filter-expression "GQ < 20"

# Select variants that PASS filters
## Site-level
### SNPs
gatk SelectVariants \
    --exclude-filtered \
    -V $output_dir/__SAMPLE___filter_snps.gatk.vcf \
    -O $output_dir/__SAMPLE___filter_site_snps.gatk.vcf
### INDELs
gatk SelectVariants \
    --exclude-filtered \
    -V $output_dir/__SAMPLE___filter_indels.gatk.vcf \
    -O $output_dir/__SAMPLE___filter_site_indels.gatk.vcf
## Sample-level
### SNPs
grep -v -E "DP_filter|GQ_filter" \
    $output_dir/__SAMPLE___filter_site_snps.gatk.vcf \
> $output_dir/__SAMPLE___filtered_snps.gatk.vcf
### INDELs
grep -v -E "DP_filter|GQ_filter" \
    $output_dir/__SAMPLE___filter_site_indels.gatk.vcf \
> $output_dir/__SAMPLE___filtered_indels.gatk.vcf

# Merge SNPs & INDELs
## Rename "SAMPLE" column in the SNP VCF file
bgzip $output_dir/__SAMPLE___filtered_snps.gatk.vcf
tabix -p vcf $output_dir/__SAMPLE___filtered_snps.gatk.vcf.gz
bcftools reheader \
    -s <(echo "__SAMPLE__ __SAMPLE___snps") \
    -o $output_dir/__SAMPLE___filtered_snps_renamed.gatk.vcf.gz \
    $output_dir/__SAMPLE___filtered_snps.gatk.vcf.gz
tabix -p vcf $output_dir/__SAMPLE___filtered_snps_renamed.gatk.vcf.gz
## Rename "SAMPLE" column in the INDEL VCF file
bgzip $output_dir/__SAMPLE___filtered_indels.gatk.vcf
tabix -p vcf $output_dir/__SAMPLE___filtered_indels.gatk.vcf.gz
bcftools reheader \
    -s <(echo "__SAMPLE__ __SAMPLE___indels") \
    -o $output_dir/__SAMPLE___filtered_indels_renamed.gatk.vcf.gz \
    $output_dir/__SAMPLE___filtered_indels.gatk.vcf.gz
tabix -p vcf $output_dir/__SAMPLE___filtered_indels_renamed.gatk.vcf.gz
## Merge SNP and INDEL VCF files
bcftools merge \
    -m all \
    -O z \
    -o $output_dir/__SAMPLE___filtered_variants.gatk.vcf.gz \
    $output_dir/__SAMPLE___filtered_snps_renamed.gatk.vcf.gz \
    $output_dir/__SAMPLE___filtered_indels_renamed.gatk.vcf.gz
tabix -p vcf $output_dir/__SAMPLE___filtered_variants.gatk.vcf.gz

# Move intermediate files to the tmp directory
mv \
$output_dir/__SAMPLE___raw_snps.gatk.vcf \
$output_dir/__SAMPLE___raw_snps.gatk.vcf.idx \
$output_dir/__SAMPLE___raw_indels.gatk.vcf \
$output_dir/__SAMPLE___raw_indels.gatk.vcf.idx \
$output_dir/__SAMPLE___filter_snps.gatk.vcf \
$output_dir/__SAMPLE___filter_snps.gatk.vcf.idx \
$output_dir/__SAMPLE___filter_indels.gatk.vcf \
$output_dir/__SAMPLE___filter_indels.gatk.vcf.idx \
$output_dir/__SAMPLE___filter_site_snps.gatk.vcf \
$output_dir/__SAMPLE___filter_site_snps.gatk.vcf.idx \
$output_dir/__SAMPLE___filter_site_indels.gatk.vcf \
$output_dir/__SAMPLE___filter_site_indels.gatk.vcf.idx \
$output_dir/__SAMPLE___filtered_snps.gatk.vcf.gz \
$output_dir/__SAMPLE___filtered_snps.gatk.vcf.gz.tbi \
$output_dir/__SAMPLE___filtered_snps_renamed.gatk.vcf.gz \
$output_dir/__SAMPLE___filtered_snps_renamed.gatk.vcf.gz.tbi \
$output_dir/__SAMPLE___filtered_indels.gatk.vcf.gz \
$output_dir/__SAMPLE___filtered_indels.gatk.vcf.gz.tbi \
$output_dir/__SAMPLE___filtered_indels_renamed.gatk.vcf.gz \
$output_dir/__SAMPLE___filtered_indels_renamed.gatk.vcf.gz.tbi \
$output_dir/__SAMPLE___filtered_variants.gatk.vcf.gz.tbi \
$output_dir/tmp

# Combine the SNP and INDEL columns in the merged vcf file
cp $output_dir/__SAMPLE___filtered_variants.gatk.vcf.gz $output_dir/tmp
gunzip $output_dir/__SAMPLE___filtered_variants.gatk.vcf.gz
cp $output_dir/__SAMPLE___filtered_variants.gatk.vcf $output_dir/__SAMPLE___filtered_variants.gatk.tmp
python combine_vcf_columns.py $output_dir/__SAMPLE___filtered_variants.gatk.tmp $output_dir/__SAMPLE___filtered_variants.gatk.vcf __SAMPLE__
rm $output_dir/__SAMPLE___filtered_variants.gatk.tmp



# ---------
# FreeBayes
# ---------

# Call variants
freebayes \
    -f $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    $output_dir/__SAMPLE___sorted_dedup_reads.bam \
> $output_dir/__SAMPLE___raw_variants.freebayes.vcf

# Filter variants
bcftools filter \
    -i 'QUAL >= 20 && FORMAT/DP >= 10' \
    $output_dir/__SAMPLE___raw_variants.freebayes.vcf \
    -o $output_dir/__SAMPLE___filtered_variants.freebayes.vcf



# ----------------
# Bcftools mpileup
# ----------------

# Call variants
## Set -d (--max-depth) to an extremely high value to "disable" the depth limit
bcftools mpileup \
    -Ou \
    -d 100000 \
    -f $output_dir/reference/${ref_fasta_basename}.${ref_fasta_extension} \
    --annotate DP,AD \
    --threads $OMP_NUM_THREADS \
    $output_dir/__SAMPLE___sorted_dedup_reads.bam \
| bcftools call \
    -mv -Ob \
    -o $output_dir/__SAMPLE___raw_variants.bcftools.bcf
## Convert bcf to vcf
bcftools view -Ov -o $output_dir/__SAMPLE___raw_variants.bcftools.vcf $output_dir/__SAMPLE___raw_variants.bcftools.bcf

# Filter variants
bcftools filter \
    -i 'QUAL >= 20 && FORMAT/DP >= 10' \
    -Ov \
    $output_dir/__SAMPLE___raw_variants.bcftools.bcf \
    -o $output_dir/__SAMPLE___filtered_variants.bcftools.vcf

# Move intermediate files to the tmp directory
mv \
$output_dir/__SAMPLE___raw_variants.bcftools.bcf \
$output_dir/tmp



######### Summary #########

# Print some messages
echo "GATK HaplotypeCaller detected $(grep -v '^#' $output_dir/__SAMPLE___filtered_variants.gatk.vcf | wc -l) variants after filtering"
echo "FreeBayes detected $(grep -v '^#' $output_dir/__SAMPLE___filtered_variants.freebayes.vcf | wc -l) variants after filtering"
echo "Bcftools mpileup detected $(grep -v '^#' $output_dir/__SAMPLE___filtered_variants.bcftools.vcf | wc -l) variants after filtering"

# Move the tmp directory to scratch
mv $output_dir/tmp $tmp_dir
echo "Intermediate files could be found in $tmp_dir"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___bwa_gatk_freebayes_bcftools-${SLURM_JOB_ID}.out
    fi
fi
