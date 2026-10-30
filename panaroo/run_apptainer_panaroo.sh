#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_panaroo

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <container> <input_dir_path> <gff_extension> <output_dir_path>"
    echo "Usage on cluster: sbatch $0 <container> <input_dir_path> <gff_extension> <output_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
container="$1"
input_dir_path="$2"
gff_extension="$3"
output_dir_path="$4"

# Check if the input directory exists
if [ ! -d "$input_dir_path" ]; then
    echo "Error: The $input_dir_path directory does not exist"
    exit 1
fi

# Check if the input directory contains gff files
if ! ls "$input_dir_path"/*.${gff_extension} 1> /dev/null 2>&1; then
    echo "Error: No $gff_extension files found in the $input_dir_path directory"
    exit 1
fi



######### Pangenome clustering #########

apptainer exec \
-W $TMPDIR \
$container panaroo \
--input $input_dir_path/*.${gff_extension} \
--out_dir $output_dir_path \
--clean-mode strict \
--alignment core \
--remove-invalid-genes \
--threads $OMP_NUM_THREADS 



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_panaroo-${SLURM_JOB_ID}.out
    fi
fi
