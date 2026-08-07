#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_pyseer_burden_lmm

######### Preprocessing #########

# Check if the correct number of arguments is provided
# Number of mandatory and optional arguments
MANDATORY_ARGS=6
TOTAL_ARGS=8

if [ $# -lt $MANDATORY_ARGS ]; then
	echo "Error: You must provide at least $MANDATORY_ARGS arguments."
    echo "Usage on terminal: bash $0 <container> <pyseer_script_dir> <vcf_file> <vcf_region> <phenotype_file> <phylogenetic_tree> [alpha] [prefix]"
    echo "Usage on cluster: sbatch $0 <container> <pyseer_script_dir> <vcf_file> <vcf_region> <phenotype_file> <phylogenetic_tree> [alpha] [prefix]"
    exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
	echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."
    echo "Usage on terminal: bash $0 <container> <pyseer_script_dir> <vcf_file> <vcf_region> <phenotype_file> <phylogenetic_tree> [alpha] [prefix]"
    echo "Usage on cluster: sbatch $0 <container> <pyseer_script_dir> <vcf_file> <vcf_region> <phenotype_file> <phylogenetic_tree> [alpha] [prefix]"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.4.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=$HOME/scratch/

# Asign arguments to variables
container="$1"
pyseer_script_dir="$2"
vcf_file="$3"
vcf_region="$4"
phenotype_file="$5"
phylogenetic_tree="$6"
alpha="${7:-0.05}"
prefix="${8:-pyseer_lmm}"

# Check vcf index file
if [ ! -f "${vcf_file}.tbi" ]; then
    echo "Error: VCF index file ${vcf_file}.tbi not found."
    exit 1
fi

# Declare variables
timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="pyseer_lmm_results_${timestamp}"
logfile="$output_dir/pyseer_run.log"

# Make output directory
if [ -d "$output_dir" ]; then
    mv $output_dir ${output_dir}_old_${timestamp}
fi
mkdir -p $output_dir

# List of expected output files
output_files=(
    "${prefix}_similarity_matrix.tsv"
    "${prefix}_all_variants.tsv"
    "${prefix}_patterns.txt"
    "${prefix}_bonferroni_threshold.txt"
    "${prefix}_significant_variants.tsv"
    "${prefix}_qq_plot.png"
)

# Write log file
exec > >(tee "$logfile") 2>&1
echo "Working directory: $PWD"
echo "Command: $0 $*"



######### Pyseer #########
# Linear mixed model (FaST-LMM)

# Calculate similarity matrix from phylogenetic tree
echo "Calculating similarity matrix from phylogenetic tree..."
apptainer run \
-W $TMPDIR \
$container python \
$pyseer_script_dir/phylogeny_distance.py --lmm \
$phylogenetic_tree \
> ${prefix}_similarity_matrix.tsv

# Run pyseer with linear mixed model
echo "Running pyseer with linear mixed model..."
apptainer run \
-W $TMPDIR \
$container pyseer --lmm \
--phenotypes $phenotype_file \
--vcf $vcf_file \
--burden $vcf_region \
--similarity ${prefix}_similarity_matrix.tsv \
--output-patterns ${prefix}_patterns.txt \
--cpu $OMP_NUM_THREADS \
> ${prefix}_all_variants.tsv
sed -i '/^$/d' ${prefix}_all_variants.tsv # remove empty lines

# qq-plot
echo "Generating qq-plot..."
apptainer run \
-W $TMPDIR \
$container python \
$pyseer_script_dir/qq_plot.py \
${prefix}_all_variants.tsv
mv qq_plot.png ${prefix}_qq_plot.png

# Calculate Bonferroni threshold
echo "Calculating Bonferroni threshold with alpha=$alpha..."
apptainer run \
-W $TMPDIR \
$container python \
$pyseer_script_dir/count_patterns.py \
--alpha $alpha \
${prefix}_patterns.txt \
> ${prefix}_bonferroni_threshold.txt
threshold=$(sed -n '2p' ${prefix}_bonferroni_threshold.txt | cut -f2) 
echo "Bonferroni threshold: $threshold"

# Get significant variants
echo "Getting significant variants..."
cat <(head -1 ${prefix}_all_variants.tsv) \
<(awk -v threshold=$threshold '$4<threshold {print $0}' ${prefix}_all_variants.tsv)\
> ${prefix}_significant_variants.tsv
echo "Pyseer with linear mixed model detected $(tail -n +2 ${prefix}_significant_variants.tsv | wc -l) significant variants"

# Sort significant variants by lrt-pvalue
echo "Sorting significant variants by lrt-pvalue..."
sort -k4 -g ${prefix}_significant_variants.tsv > ${prefix}_significant_variants_sorted.tsv
mv ${prefix}_significant_variants_sorted.tsv ${prefix}_significant_variants.tsv

# Move results to output directory
for file in "${output_files[@]}"
do
    if [ -f "$file" ]; then
        mv $file $output_dir
    fi
done
echo "Output files saved to $output_dir"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_pyseer_burden_lmm-${SLURM_JOB_ID}.out
    fi
fi
