#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <batch_size> <scheduler> <container> <fasta_dir_path> <fasta_extension>"
    exit 1
fi

# Asign arguments to variables
template="$1"
sample_list="$2"
batch_size="$3"
scheduler="$4"
container="$5"
fasta_dir_path="$6"
fasta_extension="$7"

# Set variables
subsample_list="sub"$sample_list
massiveSubmit=$template"-Launch.sh"
sample_list_length=$(wc -l < $sample_list)
nb_subsample=$((($sample_list_length / $batch_size) + 1))
suffix_length="${#nb_subsample}"

# Clear or create files
> "$subsample_list"
> "$massiveSubmit"

# Clean up old batch files if they exist
rm -rf ${sample_list}_batch_*

# Split sample list into batches
split -a $suffix_length -l $batch_size -d $sample_list ${sample_list}_batch_

# Append the names of the generated files to the subsample list
for file in ${sample_list}_batch_*; do
    echo "$file" >> "$subsample_list"
done

echo "Splitting complete! Subsamples created with prefix ${sample_list}_batch_"
echo "Subsample lists saved to $subsample_list"

# Split samples into sub directory 
for X in $(cat $subsample_list)
do
#    rm -rf ${X}_dir
    mkdir ${X}_dir
    for Y in $(cat $X)
    do
        cp ${fasta_dir_path}/${Y}.${fasta_extension} ${X}_dir    
    done
    echo "Input directory created for $X"
done

# Script cloning
for i in $(cat $subsample_list)
do
    item=$template-$i.sh
    sed "s/__SAMPLELIST__/$i/g"	$template > $item
    echo $scheduler $item $container ${i}_dir $fasta_extension ${i}_checkm2_result >> $massiveSubmit
done
