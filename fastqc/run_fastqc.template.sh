#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___fastqc

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <fastq_dir_path> <output_dir_path>"
    echo "Usage on cluster: sbatch $0 <fastq_dir_path> <output_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 fastqc/0.12.1

# Asign arguments to variables
fastq_dir_path="$1"
output_dir_path="$2"

# Declare variables
R1_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
R2_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)

# Check if the fastq file paths are correct
fastq_paths=("$R1_path" "$R2_path")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done

# Check if the output directory exists
if [ ! -d "$output_dir_path" ]; then
    echo "Error: Folder $output_dir_path does not exist"
    exit 1
fi



######### FastQC #########

fastqc \
--outdir $output_dir_path \
$R1_path \
$R2_path



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___fastqc-${SLURM_JOB_ID}.out
    fi
fi
