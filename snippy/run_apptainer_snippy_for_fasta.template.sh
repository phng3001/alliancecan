#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLELIST___snippy

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container> <reference> <fasta_dir_path>"
    echo "Usage on cluster: sbatch $0 <container> <reference> <fasta_dir_path>"
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
reference="$2"
fasta_dir_path="$3"

# Declare variables
minfrac=0.9
fasta_extension=fasta



######### SNP calling #########

for X in $(cat __SAMPLELIST__)
do
    apptainer exec \
    -W $TMPDIR \
    $container snippy \
    --reference $reference \
    --ctgs $fasta_dir_path/$X/${X}.${fasta_extension} \
    --outdir $X \
    --prefix $X \
    --minfrac $minfrac \
    --cpus $OMP_NUM_THREADS \
    --tmpdir $TMPDIR \
    --force
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLELIST___snippy-${SLURM_JOB_ID}.out
    fi
fi
