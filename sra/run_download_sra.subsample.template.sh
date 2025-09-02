#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLELIST___sra

######### Preprocessing #########

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9

# Declare variables
tmp_dir=$HOME/scratch/sra

# Check if the temporary directory exists, if not create it
if [ -d "$tmp_dir" ]; then
    echo "Temporary directory $tmp_dir exists"
else
    echo "Creating temporary directory $tmp_dir"
    mkdir -p $tmp_dir
fi



######### Download reads from SRA #########

for X in $(cat __SAMPLELIST__)
do
    if [ -d "$X" ]; then
        mv "$X" "$tmp_dir"
        echo "Pre-existed folder $X has been moved to $tmp_dir"
    fi
    mkdir $X
    cd $X
    time fastq-dump --gzip -I --split-3 $X
    cd ..
    if [ -z "$(find $X -mindepth 1 -print -quit)" ]; then
        echo "Download SRA sequences failed for $X"
    fi
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLELIST___sra-${SLURM_JOB_ID}.out
    fi
fi
