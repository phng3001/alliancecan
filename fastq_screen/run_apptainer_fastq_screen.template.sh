#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___fastq_screen

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage on terminal: bash $0 <container> <conf_file> <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <container> <conf_file> <fastq_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
container="$1"
conf_file="$2"
fastq_dir_path="$3"

# Declare variables
R1_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
# R2_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)

# Check if the fastq file paths are correct
# fastq_paths=("$R1_path" "$R2_path")
fastq_paths=("$R1_path")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done



######### FastQ-Screen #########

# R1
apptainer exec \
-W $TMPDIR \
$container \
fastq_screen \
--conf $conf_file \
--threads $OMP_NUM_THREADS \
$R1_path

# # R2
# apptainer exec \
# -W $TMPDIR \
# $container \
# fastq_screen \
# --conf $conf_file \
# --threads $OMP_NUM_THREADS \
# $R2_path



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___fastq_screen-${SLURM_JOB_ID}.out
    fi
fi
