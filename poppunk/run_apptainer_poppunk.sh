#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_poppunk

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container> <input_path_file> <output_dir>"
    echo "Usage on cluster: sbatch $0 <container> <input_path_file> <output_dir>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.4.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# PopPUNK database
poppunk_data=/project/def-mouellet/Scripts_MOU/PNP/databases/poppunk/streptococcus_pneumoniae/GPS_v11
distances=$poppunk_data/GPS_v11.dists
external_clustering=$poppunk_data/GPS_v11_external_clusters.csv

# Asign arguments to variables
container="$1"
input_path_file="$2"
output_dir="$3"

# Declare variables
tmp_dir=$HOME/scratch/poppunk/$output_dir

# Check if the temporary directory exists and reset it
if [ -d "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir

# Check if the output directory exists, if yes move it to the temporary directory
if [ -d "$output_dir" ]; then
    mv "$output_dir" "$tmp_dir"
    echo "Pre-existed folder $output_dir has been moved to $tmp_dir"
fi
mkdir $output_dir



######### PopPUNK assignment #########

apptainer exec \
-W $TMPDIR \
$container poppunk_assign \
--db $poppunk_data \
--distances $distances \
--external-clustering $external_clustering \
--query $input_path_file \
--output $output_dir \
--threads $OMP_NUM_THREADS



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_poppunk-${SLURM_JOB_ID}.out
    fi
fi
