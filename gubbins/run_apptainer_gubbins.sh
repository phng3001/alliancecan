#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_gubbins

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage on terminal: bash $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed>"
    echo "Usage on cluster: sbatch $0 <container> <core_alignment> <tree_builder {fasttree,raxml,iqtree}> <n_iteration> <prefix> <seed>" 
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



######### Recombination detection #########

apptainer exec \
-W $TMPDIR \
$container run_gubbins.py \
--tree-builder $tree_builder \
--sh-test \
--iterations $n_iteration \
--prefix $prefix \
--threads $OMP_NUM_THREADS \
--seed $seed \
$core_alignment



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_gubbins-${SLURM_JOB_ID}.out
    fi
fi
