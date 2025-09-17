#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <scheduler> <fastq_dir_path> <adapter_file_path>"
    exit 1
fi

# Parse the arguments
template="$1"
sample_list="$2"
scheduler="$3"
fastq_dir_path="$4"
adapter_file_path="$5"

# Set variables
massiveSubmit=$template"-Launch.sh"

# Clear or create files
> "$massiveSubmit"

# Script cloning
for i in $(cat $sample_list)
do
    item=$template-$i.sh
    sed "s/__SAMPLE__/$i/g" $template > $item
    echo $scheduler $item $fastq_dir_path $adapter_file_path >> $massiveSubmit
done
