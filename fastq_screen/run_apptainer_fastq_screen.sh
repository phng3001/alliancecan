#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_fastq_screen

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container> <conf_file> <fastq_file>"
    echo "Usage on cluster: sbatch $0 <container> <conf_file> <fastq_file>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
container="$1"
conf_file="$2"
fastq_file="$3"



######### FastQ-Screen #########

apptainer exec \
-W $TMPDIR \
$container \
fastq_screen \
--conf $conf_file \
--threads $OMP_NUM_THREADS \
$fastq_file



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_fastq_screen-${SLURM_JOB_ID}.out
    fi
fi
