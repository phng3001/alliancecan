#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_snippy_fasttree

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container> <recombination_free_alignment> <output_prefix>"
    echo "Usage on cluster: sbatch $0 <container> <recombination_free_alignment> <output_prefix>"
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
alignment="$2"
prefix="$3"



######### Core SNP phylogeny #########

apptainer exec \
-W $TMPDIR \
$container snp-sites \
-c $alignment \
> ${prefix}.aln

if [ -f "${prefix}.aln" ]; then
    apptainer exec \
    $container FastTree \
    -gtr -nt ${prefix}.aln \
    > ${prefix}.nwk
else
    echo "Problem identifying SNP sites from alignment, file ${prefix}.aln was not generated."
fi



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_snippy_fasttree-${SLURM_JOB_ID}.out
    fi
fi
