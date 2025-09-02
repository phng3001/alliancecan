#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLELIST___rename_sra_read

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

# Declare variables
sra_prefix=ERR
R1_input_fastq_suffix=_1.fastq
R2_input_fastq_suffix=_2.fastq
R1_output_fastq_suffix=_R1_paired.fastq
R2_output_fastq_suffix=_R2_paired.fastq



######### Rename SRA reads #########

for X in $(cat __SAMPLELIST__)
do
    if [ -d "$X" ]; then
        cd $X
        gunzip *.fastq.gz
        python ../rename_sra_read.py $sra_prefix ${X}${R1_input_fastq_suffix} ${X}${R1_output_fastq_suffix}
        python ../rename_sra_read.py $sra_prefix ${X}${R2_input_fastq_suffix} ${X}${R2_output_fastq_suffix}
        gzip *.fastq
        echo "Reads renamed for $X"
        cd ../
    else
        echo "Warning: Folder $X does not exist"
    fi
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLELIST___rename_sra_read-${SLURM_JOB_ID}.out
    fi
fi