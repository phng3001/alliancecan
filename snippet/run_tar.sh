#!/bin/bash
#SBATCH --account=def-mouellet
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_tar

# stop on errors, stop on undefined variables, stop on broken pipelines
set -euo pipefail

# enable extglob for stripping slashes
shopt -s extglob

# Usage help
usage() {
    echo "Usage: bash/sbatch $0 <dir/file_to_archive> <compression_method>"
    echo "  compression options: gz | bz2 | xz | none"
    exit 1
}

# Check arguments
if [[ $# -ne 2 ]]; then
    usage
fi

# Assign arguments to variables
target_archive=$1
compression=$2

# Remove all trailing slashes if exist
target_archive_clean=${target_archive%%+(/)}

# Ensure input exists
if [[ ! -e "$target_archive_clean" ]]; then
    echo "Error: '$target_archive_clean' does not exist"
    exit 1
fi

# Extract just the base name (no path, no trailing slash)
archive_basename=$(basename "$target_archive_clean")



# Pick compression options
case "$compression" in
    gz)   ext="tar.gz";  tar_opts="-cvzf" ;;
    bz2)  ext="tar.bz2"; tar_opts="-cvjf" ;;
    xz)   ext="tar.xz";  tar_opts="-cvJf" ;;
    none) ext="tar";     tar_opts="-cvf"  ;;
    *)    echo "Error: Unknown compression method '$compression'"; usage ;;
esac

# Define output archive name
archive="${archive_basename}.${ext}"

# Create archive
tar $tar_opts "$archive" "$target_archive_clean"

echo "Archive $archive created"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out ${archive}-${SLURM_JOB_ID}.out
    fi
fi