#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___download_sra

######### Preprocessing #########

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 sra-toolkit/3.0.9

# Asign arguments to variables
run_id="$1"

# Declare variables
tmp_dir=$HOME/scratch/sra
run_id=__SAMPLE__

# Check if the temporary directory exists, if not create it
if [ -d "$tmp_dir" ]; then
    echo "Temporary directory $tmp_dir exists"
else
    echo "Creating temporary directory $tmp_dir"
    mkdir -p $tmp_dir
fi



######### Download reads from SRA #########

if [ -d "$run_id" ]; then
    mv "$run_id" "$tmp_dir"
    echo "Pre-existed folder $run_id has been moved to $tmp_dir"
fi

mkdir $run_id

cd $run_id

time fastq-dump --gzip -I --split-3 $run_id

cd ..

if [ -z "$(find $run_id -mindepth 1 -print -quit)" ]; then
    echo "Download SRA sequences failed for $run_id"
fi



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___download_sra-${SLURM_JOB_ID}.out
    fi
fi
