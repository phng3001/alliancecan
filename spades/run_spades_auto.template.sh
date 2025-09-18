#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=__SAMPLE___spades

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage on terminal: bash $0 <fastq_dir_path>"
    echo "Usage on cluster: sbatch $0 <fastq_dir_path>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 spades/4.0.0

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Asign arguments to variables
fastq_dir_path="$1"

# Declare variables
output_dir=__SAMPLE__
R1_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R1_paired.fastq.gz)
R2_path=$(realpath $fastq_dir_path/__SAMPLE__/__SAMPLE___R2_paired.fastq.gz)
min_contig_lgth=500 # minimum contig length
memory_limit=$((OMP_NUM_THREADS * 4))
tmp_dir=$HOME/scratch/${output_dir}_spades

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

# Check if the fastq file paths are correct
fastq_paths=("$R1_path" "$R2_path")
for file in "${fastq_paths[@]}"
do
    if [ ! -f "$file" ]; then
    echo "Error: $file does not exist"
    exit 1
    fi
done

# List of expected temporary files
tmp_files=(assembly_graph_after_simplification.gfa \
assembly_graph_after_simplification.gfa \
assembly_graph.fastg \
assembly_graph_with_scaffolds.gfa \
before_rr.fasta \
contigs.paths \
corrected \
dataset.info \
input_dataset.yaml \
misc \
mismatch_corrector \
pipeline_state \
run_spades.sh \
run_spades.yaml \
scaffolds.paths \
)



######### Assembling #########

spades.py \
--careful \
--pe1-1 $R1_path \
--pe1-2 $R2_path \
-t $OMP_NUM_THREADS \
-m $memory_limit \
--tmp-dir $tmp_dir \
-o $output_dir



######### Postprocessing #########

if [ -f "$output_dir/contigs.fasta" ]; then
    # Filter out contigs < min_contig_lgth
    awk -v min_len="$min_contig_lgth" '
        /^>/ {
            if (seq) {
                if (length(seq) >= min_len)
                    print hdr "\n" seq
            }
            hdr = $0
            seq = ""
        }
        !/^>/ {
            seq = seq $0
        }
        END {
            if (seq && length(seq) >= min_len)
                print hdr "\n" seq
        }
    ' "$output_dir/contigs.fasta" > "$output_dir/contigs_${min_contig_lgth}.fasta"

    echo "Filtered FASTA written to $output_dir/contigs_${min_contig_lgth}.fasta"
    
    # Count how many contigs were kept
    num_contigs=$(grep -c '^>' "$output_dir/contigs_${min_contig_lgth}.fasta")
    echo "$num_contigs contigs ≥ $min_contig_lgth bp were kept"
else
    echo "SPAdes failed: contigs.fasta not found in $output_dir"
fi

# Move intermediate files to the temporary directory
for file in "${tmp_files[@]}"; do
    full_path="$output_dir/$file"
    if [ -e "$full_path" ]; then
        mv "$full_path" "$tmp_dir"
    fi
done

shopt -s nullglob # enable glob expansion to safely skip unmatchable patterns
for match in "$output_dir"/K*; do
    if [ -e "$match" ]; then
        mv "$match" "$tmp_dir"
    fi
done

echo "Intermediate files could be found in $tmp_dir"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out __SAMPLE___spades-${SLURM_JOB_ID}.out
    fi
fi
