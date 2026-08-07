#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_gubbins

######### Preprocessing #########

# Check if the correct number of arguments is provided
# Number of mandatory and optional arguments
MANDATORY_ARGS=7
TOTAL_ARGS=8

if [ $# -lt $MANDATORY_ARGS ]; then
	echo "Error: You must provide at least $MANDATORY_ARGS arguments."
    echo "Usage on terminal: bash $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed> <reference_name> [min_snps]"
    echo "Usage on cluster: sbatch $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed> <reference_name> [min_snps]" 
    exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
	echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."
    echo "Usage on terminal: bash $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed> <reference_name> [min_snps]"
    echo "Usage on cluster: sbatch $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed> <reference_name> [min_snps]" 
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
# export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}
# I don't know why but the above line causes error while running on compute nodes
# e.g. "--tmpdir '/localscratch/phng3001.40951974.0' is not a directory"
export TMPDIR=$HOME/scratch/

# Asign arguments to variables
container="$1"
core_alignment="$2"
tree_builder="$3"
n_iteration="$4"
prefix="$5"
seed="$6"
reference_name="$7"
min_snps="${8:-3}"

# Declare variables
timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="gubbins_results_${tree_builder}_min_snps_${min_snps}_${timestamp}"
logfile="$output_dir/gubbins_run.log"

# Make output directory
if [ -d "$output_dir" ]; then
    mv $output_dir ${output_dir}_old
fi
mkdir -p $output_dir
echo "Output directory created: $output_dir"

# List of expected output files
output_files=(
    "${prefix}.branch_base_reconstruction.embl"
    "${prefix}.filtered_polymorphic_sites.fasta"
    "${prefix}.filtered_polymorphic_sites.phylip"
    "${prefix}.final_SH_support_tree.tre"
    "${prefix}.final_tree.tre"
    "${prefix}.log"
    "${prefix}.node_labelled.final_tree.tre"
    "${prefix}.per_branch_statistics.csv"
    "${prefix}.recombination_predictions.embl"
    "${prefix}.recombination_predictions.gff"
    "${prefix}.summary_of_snp_distribution.vcf"
)

# Write log file
exec > >(tee "$logfile") 2>&1
echo "Working directory: $PWD"
echo "Command: $0 $*"



######### Gubbins #########
# Recombination detection

apptainer exec \
-W $TMPDIR \
$container run_gubbins.py \
--min-snps $min_snps \
--tree-builder $tree_builder \
--sh-test \
--iterations $n_iteration \
--prefix $prefix \
--threads $OMP_NUM_THREADS \
--seed $seed \
$core_alignment

if [ -s "${prefix}.final_tree.tre" ]; then
    sed -i "s/Reference/$reference_name/g" ${prefix}.final_tree.tre
fi

if [ -s "${prefix}.final_SH_support_tree.tre" ]; then
    sed -i "s/Reference/$reference_name/g" ${prefix}.final_SH_support_tree.tre
fi

if [ -s "${prefix}.node_labelled.final_tree.tre" ]; then
    sed -i "s/Reference/$reference_name/g" ${prefix}.node_labelled.final_tree.tre
fi

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
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_gubbins-${SLURM_JOB_ID}.out
    fi
fi
