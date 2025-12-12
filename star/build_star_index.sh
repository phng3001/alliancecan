#!/bin/bash
# P=NP
#SBATCH --account=def-bpapad
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=build_star_index

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <index_dir_name> <ref_fasta_path>"
    echo "Usage on cluster: sbatch $0 <index_dir_name> <ref_fasta_path>"
    exit 1
fi

# Asign arguments to variables
index_dir_name="$1"
ref_fasta_path="$2"

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Load modules
module purge
module load StdEnv/2023 star/2.7.11b



######### STAR indexing #########

STAR \
--runThreadN $OMP_NUM_THREADS \
--runMode genomeGenerate \
--genomeDir $index_dir_name \
--genomeFastaFiles $ref_fasta_path \
# --genomeSAindexNbases 11



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out build_star_index-${SLURM_JOB_ID}.out
    fi
fi