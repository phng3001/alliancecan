#!/bin/bash
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_seqtk_subsampling

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <fastq_dir_path> <sample_list> <fraction> <random_seed>"
    echo "Usage on cluster: sbatch $0 <fastq_dir_path> <sample_list> <fraction> <random_seed>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 seqtk/1.4

# Export variables


# Asign arguments to variables
fastq_dir_path="$1"
sample_list="$2"
fraction="$3"
random_seed="$4"

# Declare variables
fastq_dir_path=$(realpath $fastq_dir_path)
tmp_dir=$HOME/scratch/seqtk_subsampling

# Make temporary directory
## Check if the temporary directory exists and reset it
if [ -d "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir



######### Subsampling #########

for X in $(cat $sample_list)
do
    [ -d "$X" ] && mv $X $tmp_dir && echo "Pre-existing folder $X moved to $tmp_dir"
    
    mkdir $X && cd $X
    
    seqtk sample -s{random_seed} $fastq_dir_path/$X/${X}_R1_paired.fastq.gz $fraction | gzip \
    > ${X}_${fraction}_R1_paired.fastq.gz

    seqtk sample -s{random_seed} $fastq_dir_path/$X/${X}_R2_paired.fastq.gz $fraction | gzip \
    > ${X}_${fraction}_R2_paired.fastq.gz

    cd ..

    echo "$X subsampled at $fraction fraction"

done    



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_seqtk_subsampling-${SLURM_JOB_ID}.out
    fi
fi