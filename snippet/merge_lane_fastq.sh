#!/bin/bash
# P=NP
# fastq lane merging per read direction
#SBATCH --account=def-bpapad
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=merge_lane_fastq

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <sample_list> <input_dir_path>"
    echo "Usage on cluster: sbatch $0 <sample_list> <input_dir_path>"
    exit 1
fi

# Asign arguments to variables
sample_list="$1"
input_dir_path="$2"



######### Main #########

for X in $(cat $sample_list)
do
    echo "Start processing sample $X"
    # Merge R1
    cat $input_dir_path/${X}_S*_L00*_R1_*.fastq.gz > ${X}_R1.fastq.gz
    # Merge R2
    cat $input_dir_path/${X}_S*_L00*_R2_*.fastq.gz > ${X}_R2.fastq.gz

    if [ -s "${X}_R1.fastq.gz" ] && [ -s "${X}_R2.fastq.gz" ]; then
        echo "FASTQ lane merging finished for sample $X"
    else
        echo "Error while processing sample $X"
    fi
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out merge_lane_fastq-${SLURM_JOB_ID}.out
    fi
fi