#!/bin/bash
#SBATCH --account=def-mouellet
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_tar_extract

# stop on errors, stop on undefined variables, stop on broken pipelines
set -euo pipefail

# enable extglob for stripping slashes
shopt -s extglob

# Usage help
usage() {
    echo "Usage: bash/sbatch $0 <archive_to_extract> <compression_method>"
    echo "  compression options: gz | bz2 | xz | none"
    exit 1
}

# Check arguments
if [[ $# -ne 2 ]]; then
    usage
fi

# Assign arguments to variables
target_archive=$1
compression_method=$2

# Remove all trailing slashes if exist
target_archive_clean=${target_archive%%+(/)}

# Ensure input exists
if [[ ! -e "$target_archive_clean" ]]; then
    echo "Error: '$target_archive_clean' does not exist"
    exit 1
fi



# Pick compression options
case "$compression_method" in
    gz)   ext="tar.gz";  tar_opts="-xvzf" ;;
    bz2)  ext="tar.bz2"; tar_opts="-xvjf" ;;
    xz)   ext="tar.xz";  tar_opts="-xvJf" ;;
    none) ext="tar";     tar_opts="-xvf"  ;;
    *)    echo "Error: Unknown compression method '$compression_method'"; usage ;;
esac

# Extract archive
tar $tar_opts "$target_archive_clean"

echo "Archive $target_archive_clean extracted"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out ${target_archive_clean}-${SLURM_JOB_ID}.out
    fi
fi
