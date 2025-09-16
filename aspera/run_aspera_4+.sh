#!/bin/bash
# P=NP
#SBATCH --account=def-mouellet
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name=run_aspera_upload

######### Preprocessing #########

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage on terminal: bash $0 <aspera_key_path> <upload_dir_path>"
    echo "Usage on cluster: sbatch $0 <aspera_key_path> <upload_dir_path>"
    exit 1
fi

# Asign arguments to variables
aspera_key="$1"
upload_dir="$2"

# Check if the upload directory exists
if [ ! -d "$upload_dir" ]; then
    echo "Error: Folder $upload_dir does not exist"
    exit 1
fi

# Add the ascp binaries to PATH
export PATH=~/.aspera/connect/bin:$PATH



######### Upload #########

ascp \
-v \
-i $aspera_key \
-QT -l100m -k1 \
-d $upload_dir \
subasp@upload.ncbi.nlm.nih.gov:uploads/nguyen-phuong.pham_crchudequebec.ulaval.ca_Ii0VBDSL



# Save stdout
if [[ -n "$SLURM_JOB_ID" && "$SLURM_JOB_ID" -ne 0 ]]; then
    sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8
    #sacct -j $SLURM_JOB_ID --format=JobID%16,Submit,Start,Elapsed,NCPUS,ExitCode,NodeList%8,MaxRSS
    if [[ -f "slurm-${SLURM_JOB_ID}.out" ]]; then
        mv slurm-${SLURM_JOB_ID}.out run_aspera_upload-${SLURM_JOB_ID}.out
    fi
fi
