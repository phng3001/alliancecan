#!/bin/bash
# P=NP
#SBATCH --account=def-bpapad
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___star

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <index_dir_path> <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <index_dir_path> <fastq_dir_path>"
    exit 1
fi

# Asign arguments to variables
index_dir_path="$1"
fastq_dir_path="$2"

# Declare variables
output_dir=__SAMPLE__
prefix=__SAMPLE__
index_dir_realpath=$(realpath $index_dir_path)
R1_realpath=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
R2_realpath=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)
tmp_dir=$HOME/scratch/${output_dir}_star

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Load modules
module purge
module load StdEnv/2023 star/2.7.11b samtools/1.22.1

# Make temporary directory
## Check if the temporary directory exists and reset it
if [ -e "$tmp_dir" ]; then
    rm -rf "$tmp_dir"
    echo "Pre-existed file/folder $tmp_dir has been reset"
fi
mkdir -p $tmp_dir

# Make output directory
## Check if the output directory exists, if yes move it to the temporary directory
if [ -e "$output_dir" ]; then
    mv "$output_dir" "$tmp_dir"
    echo "Pre-existed file/folder $output_dir has been moved to $tmp_dir"
fi
mkdir $output_dir



######### STAR #########

cd $output_dir

# Mapping
# --alignIntronMax 1: disable spliced alignment
# --outFilterMismatchNmax: mismatch limit
# --outFilterMultimapNmax: multimapping limit, if a read maps to more locations than threshold -> unmapped
# --limitBAMsortRAM: maximum available RAM (bytes) for sorting BAM. 4000000000 -> 4Gb

STAR \
--runThreadN $OMP_NUM_THREADS \
--genomeDir $index_dir_realpath \
--readFilesIn $R1_realpath $R2_realpath \
--readFilesCommand zcat \
--outFileNamePrefix $prefix \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--alignIntronMax 1 \
--outFilterMismatchNmax 3 \
--outFilterMultimapNmax 30 \
--limitBAMsortRAM 4000000000 \

# Indexing
samtools index \
${prefix}Aligned.sortedByCoord.out.bam

cd ..



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___star-${SLURM_JOB_ID}.out
    fi
fi
