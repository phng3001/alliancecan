#!/bin/bash
# P=NP
# blastn: 1 query (gene) vs a list of subjects (genomes), for each subject only keep the best hit 

# Check for the correct number of arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: bash $0 <query> <subject_dir> <subject_list> <output_prefix>"
    exit 1
fi

# Asign arguments to variables
query=$1
subject_dir=$2
subject_list=$3
prefix=$4

# Load modules
module purge
module load StdEnv/2023 gcc/12.3 blast+/2.14.1 samtools/1.22.1

# Clear or create output files
> "${prefix}_best_hits_nu.fasta"
> "${prefix}_best_hits_nu.tsv"

# Add header to the output tsv file
echo -e "qseqid\tsseqid\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsstrand\tlength\tqlen\tslen\tpident" > "${prefix}_best_hits_nu.tsv"

# Loop through each subject file
for X in $(cat $subject_list)
do
    subject_path="$subject_dir/$X"
    if [ ! -f "$subject_path" ]; then
        echo "Subject file $X does not exist. Skipping..."
        continue
    fi
    echo "Performing BLASTN of $query vs $X"
    blastn \
    -query $query \
    -subject $subject_path \
    -out ${query}_blastn_${X}.tmp \
    -outfmt "6 qseqid sseqid qstart qend sstart send evalue bitscore sstrand length qlen slen pident" \
    -perc_identity 75 \
    -qcov_hsp_perc 75 \
    -max_hsps 1 \
    -max_target_seqs 1 # best hit

    # extract fields
    sseqid=$(cut -f2 ${query}_blastn_${X}.tmp)
    sstart=$(cut -f5 ${query}_blastn_${X}.tmp)
    send=$(cut -f6 ${query}_blastn_${X}.tmp)

    if [ -n "$sseqid" ]; then
        # write hit to the output fasta file
        if [ "$sstart" -lt "$send" ]; then
            samtools faidx "$subject_path" "$sseqid:$sstart-$send" >> "${prefix}_best_hits_nu.fasta"
        else
            samtools faidx "$subject_path" "$sseqid:$send-$sstart" -i >> "${prefix}_best_hits_nu.fasta"          
        fi

        # append hit to the output tsv file
        cat ${query}_blastn_${X}.tmp >> "${prefix}_best_hits_nu.tsv"

    else
        echo "No hit found in $X"
    fi

    # remove temporary files
    rm ${query}_blastn_${X}.tmp
done


