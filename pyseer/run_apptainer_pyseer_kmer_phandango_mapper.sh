#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_pyseer_kmer_phandango_mapper

######### Preprocessing #########

# Check if the correct number of arguments is provided

if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <container> <significant_kmers> <reference_fasta> <output_file>"
    echo "Usage on cluster: sbatch $0 <container> <significant_kmers> <reference_fasta> <output_file>"
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
significant_kmers="$2"
reference_fasta="$3"
output_file="$4"

# Declare variables
timestamp=$(date +"%Y%m%d_%H%M%S")
logfile="pyseer_kmer_phandango_mapper_${timestamp}.log"

# Write log file
exec > >(tee "$logfile") 2>&1
echo "Working directory: $PWD"
echo "Command: $0 $*"



######### Pyseer phandango_mapper #########

echo "Mapping significant kmers to $reference_fasta..."
apptainer run \
-W $TMPDIR \
$container phandango_mapper \
$significant_kmers \
$reference_fasta \
$output_file

if [ -f "$output_file" ]; then
    echo "Mapping completed. Output saved to $output_file"
else
    echo "Problem mapping significant kmers, output file $output_file was not generated."
fi



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_pyseer_kmer_phandango_mapper-${SLURM_JOB_ID}.out
    fi
fi
