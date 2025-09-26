#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_prokka

######### Preprocessing #########

# Check if the number of arguments is correct
# Number of mandatory and optional arguments
MANDATORY_ARGS=2
TOTAL_ARGS=6

if [ $# -lt $MANDATORY_ARGS ]; then
	echo "Error: You must provide at least $MANDATORY_ARGS arguments."
	echo "Usage on terminal: bash $0 <sample_name> <fasta_file> [kingdom] [genus] [species] [strain]"
	echo "Usage on cluster: sbatch $0 <sample_name> <fasta_file> [kingdom] [genus] [species] [strain]"
	exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
	echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."
	echo "Usage on terminal: bash $0 <sample_name> <fasta_file> [kingdom] [genus] [species] [strain]"
	echo "Usage on cluster: sbatch $0 <sample_name> <fasta_file> [kingdom] [genus] [species] [strain]"
	exit 1
fi

# Load modules
module purge
module load StdEnv/2020 gcc/9.3.0 prokka/1.14.5 barrnap/0.9 rnammer/1.2

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Asign arguments to variables
sample_name="$1"
fasta_file="$2"
kingdom="${3:-Bacteria}"
genus="${4:-Genus}"
species="${5:-species}"
strain="${6:-strain}"

# Check if the user input matches any allowed option
kingdom_options=("Archaea" "Bacteria" "Mitochondria" "Viruses")
if [[ ! " ${kingdom_options[@]} " =~ " ${kingdom} " ]]; then
    echo "Error: '$kingdom' is invalid. Valid options are: ${kingdom_options[@]}"
    exit 1
fi

# Check if FASTA file exists and is readable
if [[ ! -s "$fasta_file" ]]; then
    echo "Error: FASTA file '$fasta_file' is missing or empty or unreadable"
    exit 1
fi



######### Annotation #########

# Debugging: Print variables before running Prokka
echo "Running Prokka with:"
echo "Sample: $sample_name"
echo "FASTA file: $fasta_file"
echo "Kingdom: $kingdom"
echo "Genus: $genus"
echo "Species: $species"
echo "Strain: $strain"
echo "Threads: $OMP_NUM_THREADS"

echo "The sample name $sample_name will be used as the output directory's name, output file name and locus tag prefixes"

# Run Prokka
prokka \
--rfam \
--compliant \
--outdir $sample_name \
--prefix $sample_name \
--locustag $sample_name \
--cpus $OMP_NUM_THREADS \
--kingdom $kingdom \
--genus $genus \
--species $species \
--strain $strain \
"$fasta_file"



# Save SLURM stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
	sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
	#sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS

	if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_checkm-${SLURM_JOB_ID}.out
    fi
fi
