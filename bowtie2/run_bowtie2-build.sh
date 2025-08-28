#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_bowtie2-build

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <ref_fasta_path> <bt2_index_basename>"
    echo "Usage on cluster: sbatch $0 <ref_fasta_path> <bt2_index_basename>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 bowtie2/2.5.4

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
# export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
ref_fasta_path="$1"
bt2_index_basename="$2"

# Declare variables
ref_fasta=$(realpath $ref_fasta_path)



######### Bowtie2 indexing #########

bowtie2-build \
$ref_fasta \
$bt2_index_basename \
--threads $OMP_NUM_THREADS



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_bowtie2-build-${SLURM_JOB_ID}.out
    fi
fi
