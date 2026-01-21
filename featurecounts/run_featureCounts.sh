#!/bin/bash
# P=NP
#SBATCH --account=def-bpapad
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_featureCounts

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: bash $0 <annotation_gtf/gff> <attribute_type> <bam_file_list> <output_tsv>"
    exit 1
fi

# Asign arguments to variables
annotation_file="$1"
attribute_type="$2"
bam_file_list="$3"
output_file="$4"

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Load modules
module purge
module load StdEnv/2023 subread/2.0.6



######### featureCounts #########

# -p paired-end
# -s strandness: 0 (unstranded, default) | 1 (stranded/forward) | 2 (stranded/reverse)
# -t feature type: e.g. transcript, exon
# -g attribute type: e.g. gene_id (default), transcript_id
# -B: Require both ends mapped
# -C: Exclude chimeric fragments
# -M --fraction: multi-mapping reads included, fractional counts assigned
# --primary: Count primary alignments

featureCounts \
-T $OMP_NUM_THREADS \
-p \
-s 2 \
-B -C \
--primary \
-g $attribute_type \
-a $annotation_file \
-o $output_file \
$(cat $bam_file_list)



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_featureCounts-${SLURM_JOB_ID}.out
    fi
fi
