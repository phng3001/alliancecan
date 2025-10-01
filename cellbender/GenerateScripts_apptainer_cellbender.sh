#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 6 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <scheduler> <container_path> <input_dir_path> <output_dir_path>"
    exit 1
fi

# Parse the arguments
template="$1"
sample_list="$2"
scheduler="$3"
container_path="$4"
input_dir_path="$5"
output_dir_path="$6"

# Set variables
massiveSubmit=$template"-Launch.sh"

# Clear or create files
> "$massiveSubmit"

# Script cloning
for i in $(cat $sample_list)
do
    item=$template-$i.sh
    sed "s/__SAMPLE__/$i/g" $template > $item
    echo $scheduler $item $container_path $input_dir_path $output_dir_path >> $massiveSubmit
done
