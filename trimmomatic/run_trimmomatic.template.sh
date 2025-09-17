#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___trimmomatic

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <fastq_dir_path> <adapter_file_path>"
    echo "Usage on cluster: sbatch $0 <fastq_dir_path> <adapter_file_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 trimmomatic/0.39

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Asign arguments to variables
fastq_dir_path="$1"
adapter_file_path="$2"

# Declare variables
output_dir=__SAMPLE__
R1_path=$(realpath $fastq_dir_path/__SAMPLE___1.fq.gz)
R2_path=$(realpath $fastq_dir_path/__SAMPLE___2.fq.gz)
adapter_file=$(realpath $adapter_file_path)
tmp_dir=$HOME/scratch/${output_dir}_trimmomatic

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

# Check if the file paths are correct
file_paths=("$R1_path" "$R2_path" "$adapter_file")
for file in "${file_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done



######### Trimming #########

# Trimmomatic
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar \
PE -phred33 \
-threads $OMP_NUM_THREADS \
$R1_path $R2_path \
$output_dir/__SAMPLE___R1_paired.fastq $output_dir/__SAMPLE___R1_unpaired.fastq \
$output_dir/__SAMPLE___R2_paired.fastq $output_dir/__SAMPLE___R2_unpaired.fastq \
ILLUMINACLIP:$adapter_file:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Compression
gzip $output_dir/__SAMPLE___R1_paired.fastq
gzip $output_dir/__SAMPLE___R1_unpaired.fastq
gzip $output_dir/__SAMPLE___R2_paired.fastq
gzip $output_dir/__SAMPLE___R2_unpaired.fastq



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___trimmomatic-${SLURM_JOB_ID}.out
    fi
fi


