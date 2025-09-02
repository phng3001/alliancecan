#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 4 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <batch_size> <scheduler>"
    exit 1
fi

# Asign arguments to variables
template="$1"
sample_list="$2"
batch_size="$3"
scheduler="$4"

# Set variables
subsample_list="sub"$sample_list
massiveSubmit=$template"-Launch.sh"

# Clear or create files
> "$subsample_list"
> "$massiveSubmit"

# Split sample list into batches
split -l $batch_size -d $sample_list ${sample_list}_batch_

# Append the names of the generated files to the subsample list
for file in ${sample_list}_batch_*; do
    echo "$file" >> "$subsample_list"
done

echo "Splitting complete! Subsamples created with prefix ${sample_list}_batch_"
echo "Subsample lists saved to $subsample_list"

# Script cloning
for i in $(cat $subsample_list)
do
    item=$template-$i.sh
    sed "s/__SAMPLELIST__/$i/g"	$template > $item
    echo $scheduler $item >> $massiveSubmit
done
