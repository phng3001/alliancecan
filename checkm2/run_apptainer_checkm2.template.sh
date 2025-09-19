#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLELIST___checkm2

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <container> <fasta_dir_path> <fasta_extension> <output_dir_name>"
    echo "Usage on cluster: sbatch $0 <container> <fasta_dir_path> <fasta_extension> <output_dir_name>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Reference data must be downloaded and decompressed in this directory
export CHECKM2DB=/project/def-mouellet/Scripts_MOU/PNP/databases/checkm2/CheckM2_database

# Asign arguments to variables
container="$1"
fasta_dir_path="$2"
fasta_extension="$3"
output_dir_name="$4"



######### CheckM2 #########

# Run checkm2 predict
apptainer exec \
-B $CHECKM2DB \
-W $TMPDIR \
$container checkm2 predict \
-x $fasta_extension \
--input $fasta_dir_path \
--output-directory $output_dir_name \
--database_path $CHECKM2DB/uniref100.KO.1.dmnd \
--threads $OMP_NUM_THREADS \
--force



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLELIST___checkm2-${SLURM_JOB_ID}.out
    fi
fi


