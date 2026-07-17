#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_orthofinder

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <container> <input_dir_path>"
    echo "Usage on cluster: sbatch $0 <container> <input_dir_path>"
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
input_dir_path="$2"



######### OrthoFinder #########

apptainer exec \
-W $TMPDIR \
$container \
orthofinder \
-f $input_dir_path \
-t $OMP_NUM_THREADS



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_orthofinder-${SLURM_JOB_ID}.out
    fi
fi
