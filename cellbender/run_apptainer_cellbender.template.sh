#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___cellbender_remove-background

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container_path> <input_dir_path> <output_dir_path>"
    echo "Usage on cluster: sbatch $0 <container_path> <input_dir_path> <output_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
container_path="$1"
input_dir_path="$2"
output_dir_path="$3"

# Declare variables
container=$(realpath $container_path)
input=$(realpath $input_dir_path/__SAMPLE__.h5ad)
output=$(realpath $output_dir_path/__SAMPLE___denoised)
tmp_dir=__SAMPLE___tmp
droplets=50000

# Make temporary directory
if [ -d "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir



######### CellBender remove-background #########

cd $tmp_dir

apptainer run --nv \
-W $TMPDIR \
$container \
cellbender remove-background \
--input $input \
--output $output \
--total-droplets-included $droplets \
--cuda

cd ..



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___cellbender_remove-background-${SLURM_JOB_ID}.out
    fi
fi
