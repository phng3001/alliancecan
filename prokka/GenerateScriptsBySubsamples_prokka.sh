#!/bin/bash
# P=NP

# Check if the number of arguments is correct
# Number of mandatory and optional arguments
MANDATORY_ARGS=6
TOTAL_ARGS=9

if [ $# -lt $MANDATORY_ARGS ]; then
	echo "Error: You must provide at least $MANDATORY_ARGS arguments."
    echo "Usage: bash $0 <template.sh> <sample_list> <batch_size> <scheduler> <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	exit 1
fi

if [ $# -gt $TOTAL_ARGS ]; then
	echo "Error: Too many arguments. You can provide a maximum of $TOTAL_ARGS arguments."
	echo "Usage: bash $0 <template.sh> <sample_list> <batch_size> <scheduler> <fasta_dir_path> <fasta_extension> [kingdom] [genus] [species]"
	exit 1
fi

# Asign arguments to variables
template="$1"
sample_list="$2"
batch_size="$3"
scheduler="$4"
fasta_dir_path="$5"
fasta_extension="$6"
kingdom="$7"
genus="$8"
species="$9"

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
	echo $scheduler $item $fasta_dir_path $fasta_extension $kingdom $genus $species >> $massiveSubmit
done