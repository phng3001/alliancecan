#!/bin/bash
#SBATCH --account=def-mouellet
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_apptainer_ncbi-genome-download_refseq_bacteria

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage on terminal: bash $0 <container> <species_taxid> <format> <output_prefix>"
    echo "Usage on cluster: sbatch $0 <container> <species_taxid> <format> <output_prefix>"
    exit 1
fi

# Load modules
module purge
module load StdEnv/2023 apptainer/1.3.5

# Export variables
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}

# Asign arguments to variables
container="$1"
species_taxid="$2"
format="$3"
output_prefix="$4"



######### Genome downloading #########

# List assemblies for download
apptainer exec \
--bind /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:/etc/pki/tls/certs/ca-bundle.crt \
$container ncbi-genome-download bacteria \
--species-taxids $species_taxid \
--dry-run > ${output_prefix}_list.tsv
# Remove the first line (message)
sed -i '1d' ${output_prefix}_list.tsv 

# Download assemblies
apptainer exec \
--bind /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem:/etc/pki/tls/certs/ca-bundle.crt \
$container ncbi-genome-download bacteria \
--species-taxids $species_taxid \
-o ${output_prefix}_genomes \
--formats $format \
--flat-output \
--metadata-table ${output_prefix}_metadata.tsv \
-p $OMP_NUM_THREADS \
--progress-bar
# Unzip files
gunzip ${output_prefix}_genomes/*.gz



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS    
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_apptainer_ncbi-genome-download_refseq_bacteria-${SLURM_JOB_ID}.out
    fi
fi
