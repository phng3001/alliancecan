#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_snippy

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage on terminal: bash $0 <container> <reference> <sample_name> <R1_path> <R2_path> <minfrac>"
    echo "Usage on cluster: sbatch $0 <container> <reference> <sample_name> <R1_path> <R2_path> <minfrac>"
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
sample_name="$3"
R1_path="$4"
R2_path="$5"
minfrac="$6"

# Declare variables
reference_basename="${reference%.*}"



######### SNP calling #########

apptainer exec \
-W $TMPDIR \
$container snippy \
--reference $reference \
--R1 $R1_path \
--R2 $R2_path \
--outdir ${sample_name}_mapping_${reference_basename} \
--prefix $sample_name \
--minfrac $minfrac \
--cpus $OMP_NUM_THREADS \
--tmpdir $TMPDIR \
--force



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_snippy_${sample_name}_mapping_${reference_basename}-${SLURM_JOB_ID}.out
    fi
fi
