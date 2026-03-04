#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_rseqc_infer_experiment

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <container> <sample_list> <bam_dir_path> <ref_bed_path>"
    echo "Usage on cluster: sbatch $0 <container> <sample_list> <bam_dir_path> <ref_bed_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}

# Asign arguments to variables
container="$1"
sample_list="$2"
bam_dir_path="$3"
ref_bed_path="$4"



######### RSeQC #########

for sample in $(cat $sample_list)
do
    bam_path=$(realpath $bam_dir_path/$sample/${sample}.bam)
    apptainer exec \
    -W $TMPDIR \
    $container \
    infer_experiment.py \
    -i $bam_path \
    -r $ref_bed_path \
    > ${sample}_rseqc.txt
done



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_rseqc_infer_experiment-${SLURM_JOB_ID}.out
    fi
fi
