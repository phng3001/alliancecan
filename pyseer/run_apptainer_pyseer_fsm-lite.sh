#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=256G
#SBATCH --job-name=run_apptainer_pyseer_fsm-lite

######### Preprocessing #########

# Check if the correct number of arguments is provided

if [ "$#" -ne 7 ]; then
    echo "Usage on terminal: bash $0 <container> <input_file_list> <output_file> <min_length> <max_length> <min_supp> <max_supp>"
    echo "Usage on cluster: sbatch $0 <container> <input_file_list> <output_file> <min_length> <max_length> <min_supp> <max_supp>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.4.5

# Export variables
export TMPDIR=$HOME/scratch/

# Asign arguments to variables
container="$1"
input_file_list="$2"
output_file="$3"
min_length="$4"
max_length="$5"
min_supp="$6"
max_supp="$7"

# Declare variables
timestamp=$(date +"%Y%m%d_%H%M%S")
tmp_index="fsm_kmers_${timestamp}"
logfile="fsm-lite_run_${timestamp}.log"

# Write log file
exec > >(tee "$logfile") 2>&1
echo "Working directory: $PWD"
echo "Command: $0 $*"



######### fsm-lite #########
# Count k-mers in the input files

# Run fsm-lite to count k-mers
echo "Running fsm-lite to count k-mers..."
apptainer run \
-W $TMPDIR \
$container fsm-lite \
--list $input_file_list \
--min $min_length \
--max $max_length \
--minsupp $min_supp \
--maxsupp $max_supp \
--tmp $tmp_index \
--verbose \
> $output_file
gzip -k $output_file

echo "fsm-lite run completed. Number of k-mers detected: $(cat $output_file | wc -l)"
echo "Output file: $output_file"
echo "Compressed output file: ${output_file}.gz"



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_pyseer_fsm-lite-${SLURM_JOB_ID}.out
    fi
fi
