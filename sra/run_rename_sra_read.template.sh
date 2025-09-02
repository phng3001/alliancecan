#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___rename_sra_read

######### Preprocessing #########

# Prerequisites
scripts=("rename_sra_read.py")

# Check if the prerequisite scripts exist in the current directory
for script in "${scripts[@]}"
do
    if [ ! -f "$script" ]; then
    echo "Error: $script not found in the current directory"
    exit 1
    fi
done

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <sra_prefix> <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <sra_prefix> <fastq_dir_path>"
    exit 1
fi

# Asign arguments to variables
sra_prefix="$1"
fastq_dir_path="$2"

# Declare variables
R1_input_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___1.fastq)
R2_input_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___2.fastq)
R1_output_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq)
R2_output_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq)

# Check if the fastq input file paths are correct
fastq_paths=("${R1_input_path}.gz" "${R2_input_path}.gz")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done



######### Rename SRA reads #########

gunzip ${R1_input_path}.gz
gunzip ${R2_input_path}.gz
python rename_sra_read.py $sra_prefix $R1_input_path $R1_output_path
python rename_sra_read.py $sra_prefix $R2_input_path $R2_output_path
gzip $R1_input_path
gzip $R2_input_path
gzip $R1_output_path
gzip $R2_output_path
echo "Renamed reads written at ${R1_output_path}.gz and ${R2_output_path}.gz"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___rename_sra_read-${SLURM_JOB_ID}.out
    fi
fi
