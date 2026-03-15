#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_eggnog_protein

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <eggnog_container> <fasta_file>"
    echo "Usage on cluster: sbatch $0 <eggnog_container> <fasta_file>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.4.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export TMPDIR=${SLURM_TMPDIR:-$HOME/scratch/}
# eggNOG-mapper databases
export EGGNOG_DATA_DIR=/project/def-mouellet/Scripts_MOU/PNP/databases/eggnog

# Asign arguments to variables
eggnog_container="$1"
fasta_file="$2"

# Declare variables
filename="${fasta_file##*/}"
output_dir="${filename%.*}"
tmp_dir=$HOME/scratch/${output_dir}_eggnog

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



######### Annotation #########

apptainer exec \
-W $TMPDIR \
$eggnog_container emapper.py \
-i $fasta_file \
-o $output_dir \
--output_dir $output_dir \
--itype proteins \
-m diamond \
--pident 40 \
--query_cover 20 \
--subject_cover 20 \
--evalue 0.001 \
--score 60 \
--dmnd_ignore_warnings \
--tax_scope auto \
--target_orthologs all \
--report_orthologs \
--go_evidence non-electronic \
--pfam_realign none \
--decorate_gff yes \
--excel \
--cpu $OMP_NUM_THREADS \
--scratch_dir $tmp_dir \
--temp_dir $tmp_dir 



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_eggnog_protein-${SLURM_JOB_ID}.out
    fi
fi