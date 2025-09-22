#!/bin/bash
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_tar

# Check for the correct number of arguments
if [ "$#" -ne 1 ]; then
    echo "Usage on terminal: bash $0 <dir_to_archive>"
    echo "Usage on terminal: sbatch $0 <dir_to_archive>"
    exit 1
fi

# Assign arguments to variables
dir_to_archive=$1



# Run tar
tar -cvzf ${dir_to_archive}.tar.gz $dir_to_archive



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out ${dir_to_archive}_tar-${SLURM_JOB_ID}.out
    fi
fi