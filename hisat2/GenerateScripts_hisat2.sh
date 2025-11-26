#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <scheduler> <ref_fasta_path> <fastq_dir_path>"
    exit 1
fi

# Parse the arguments
template="$1"
sample_list="$2"
scheduler="$3"
ref_fasta_path="$4"
fastq_dir_path="$5"

# Set variables
massiveSubmit=$template"-Launch.sh"

# Clear or create files
> "$massiveSubmit"

# Script cloning
for i in $(cat $sample_list)
do
    item=$template-$i.sh
    sed "s/__SAMPLE__/$i/g" $template > $item
    echo $scheduler $item $ref_fasta_path $fastq_dir_path >> $massiveSubmit
done
