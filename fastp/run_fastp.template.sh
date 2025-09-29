#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___fastp

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage on terminal: bash $0 <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <fastq_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 fastp/0.24.0

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Asign arguments to variables
fastq_dir_path="$1"

# Declare variables
output_dir=__SAMPLE__
R1_path=$(realpath $fastq_dir_path/__SAMPLE___R1.fastq.gz)
R2_path=$(realpath $fastq_dir_path/__SAMPLE___R2.fastq.gz)
tmp_dir=$HOME/scratch/${output_dir}_fastp

# Check if the fastq file paths are correct
fastq_paths=("$R1_path" "$R2_path")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done

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



######### Trimming #########

fastp \
-i $R1_path \
-I $R2_path \
-o $output_dir/__SAMPLE___R1_paired.fastq.gz \
-O $output_dir/__SAMPLE___R2_paired.fastq.gz \
--unpaired1 $output_dir/__SAMPLE___R1_unpaired.fastq.gz \
--unpaired2 $output_dir/__SAMPLE___R2_unpaired.fastq.gz \
--detect_adapter_for_pe \
-l 20 \
-q 20 \
-g \
--correction \
-w $OMP_NUM_THREADS \
-h $output_dir/__SAMPLE___fastp.html \
-j $output_dir/__SAMPLE___fastp.json



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___fastp-${SLURM_JOB_ID}.out
    fi
fi