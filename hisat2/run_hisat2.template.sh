#!/bin/bash
# P=NP
#SBATCH --account=def-bpapad
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___hisat2

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <ref_fasta_path> <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <ref_fasta_path> <fastq_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 hisat2/2.2.1 samtools/1.20

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
ref_fasta_path="$1"
fastq_dir_path="$2"

# Declare variables
ref_fasta_realpath=$(realpath $ref_fasta_path)
output_dir=__SAMPLE__
R1_realpath=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
R2_realpath=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)
tmp_dir=$HOME/scratch/${output_dir}_hisat2

# Make temporary directory
## Check if the temporary directory exists and reset it
if [ -d "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir

# Make output directory
## Check if the output directory exists, if yes move it to the temporary directory
if [ -e "$output_dir" ]; then
    mv "$output_dir" "$tmp_dir"
    echo "Pre-existed file/folder $output_dir has been moved to $tmp_dir"
fi
mkdir $output_dir

# Check if the file paths are correct
paths=("$R1_realpath" "$R2_realpath" "$ref_fasta_realpath")
for file in "${paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done



######### HISAT2 #########

cd $output_dir

# Index reference
## Get reference fasta file
cp $ref_fasta_realpath .
## Get ref fasta file extension
ref_fasta_extension="${ref_fasta_path##*.}"
## Get ref fasta file base name without extension
ref_fasta_basename="${ref_fasta_path##*/}"
ref_fasta_basename="${ref_fasta_basename%.*}"
## Build indexes
hisat2-build \
-p $OMP_NUM_THREADS \
${ref_fasta_basename}.${ref_fasta_extension} \
$ref_fasta_basename

# Mapping
hisat2 \
-p $OMP_NUM_THREADS \
-x $ref_fasta_basename \
--no-spliced-alignment \
-1 $R1_realpath \
-2 $R2_realpath \
-S ${output_dir}.sam \
--summary-file ${output_dir}.log

# Sorting
samtools sort \
--threads $OMP_NUM_THREADS \
-o ${output_dir}.bam \
${output_dir}.sam

# Indexing
samtools index \
${output_dir}.bam

cd ..



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___hisat2-${SLURM_JOB_ID}.out
    fi
fi
