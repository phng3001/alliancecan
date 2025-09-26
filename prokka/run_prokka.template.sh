#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLELIST___prokka

######### Preprocessing #########

# Check if the number of arguments is correct
# Number of mandatory and optional arguments
MANDATORY_ARGS=2
TOTAL_ARGS=5

if [ $# -lt $MANDATORY_ARGS ]; then
	echo "Error: You must provide at least $MANDATORY_ARGS arguments."
	echo "Usage on terminal: bash $0 <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	echo "Usage on cluster: sbatch $0 <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
	echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."
	echo "Usage on terminal: bash $0 <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	echo "Usage on cluster: sbatch $0 <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	exit 1
fi

# Load modules
module purge
module load StdEnv/2020 gcc/9.3.0 prokka/1.14.5 barrnap/0.9 rnammer/1.2

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Asign arguments to variables
fasta_dir_path="$1"
fasta_extension="$2"
kingdom="${3:-Bacteria}"
genus="${4:-Genus}"
species="${5:-species}"



######### Annotation #########

for X in $(cat __SAMPLELIST__)
do
    prokka \
    --rfam \
    --compliant \
    --outdir $X \
    --prefix $X \
    --strain $X \
    --locustag $X \
    --kingdom $kingdom \
    --genus $genus \
    --species $species \
    --cpus $OMP_NUM_THREADS \
    --force \
    ${fasta_dir_path}/${X}.${fasta_extension}
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLELIST___prokka-${SLURM_JOB_ID}.out
    fi
fi
