#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___stringtie

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <ref_gff_path> <mapping_dir_path>"
    echo "Usage on cluster: sbatch $0 <ref_gff_path> <mapping_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 stringtie/3.0.1

# Asign arguments to variables
ref_gff_path="$1"
mapping_dir_path="$2"

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# Set variables
output=__SAMPLE__
ref_gff_realpath=$(realpath $ref_gff_path)
bam_file_realpath=$(realpath $mapping_dir_path/__SAMPLE__/__SAMPLE__.bam)



######### StringTie #########

stringtie \
-p $OMP_NUM_THREADS \
-G $ref_gff_realpath \
-l $output \
-o ${output}.gtf \
$bam_file_realpath



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___stringtie-${SLURM_JOB_ID}.out
    fi
fi
