#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G

# Load modules
module load StdEnv/2023 nextflow/26.04.4
module load StdEnv/2023 apptainer/1.3.5

# Environment variables
export NXF_OPTS="-Xmx2g"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "Start running pipeline using $OMP_NUM_THREADS CPUs"

# Main
nextflow run fastp.nf --ncpus $OMP_NUM_THREADS -resume -with-report -with-timeline
# nextflow run fastp.nf --task.cpus $OMP_NUM_THREADS -resume

# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___fastp-${SLURM_JOB_ID}.out
    fi
fi
