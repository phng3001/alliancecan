#!/bin/bash
# P=NP

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: bash $0 <template.sh> <sample_list> <scheduler>"
    exit 1
fi

# Asign arguments to variables
template="$1"
sample_list="$2"
scheduler="$3"

# Set variables
massiveSubmit=$template"-Launch.sh"

# Clear or create files
> "$massiveSubmit"

# Script cloning
for i in $(cat $sample_list)
do
	item=$template-$i.sh
	sed "s/__SAMPLE__/$i/g"	$template > $item
	echo $scheduler $item >> $massiveSubmit
done
