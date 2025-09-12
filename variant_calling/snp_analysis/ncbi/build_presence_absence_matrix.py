# P=NP
import sys
import argparse
from collections import defaultdict

def parse_input_file(input_file, key_columns, sample_column):
    """Parse input file to build a dictionary of variants and their presence in samples, using user-defined key columns."""
    variants = defaultdict(lambda: defaultdict(set))
    samples = set()

    with open(input_file, 'r') as file:
        header = file.readline().strip().split('\t')
        
        # Indices for relevant columns
        column_indices = {col: header.index(col) for col in key_columns if col in header}
        sample_index = header.index(sample_column)

        for line in file:
            columns = line.strip().split('\t')
            sample = columns[sample_index]
            
            # Build a tuple with user-selected key columns
            variant_key = tuple(columns[column_indices[col]] for col in key_columns)
            
            # Add this variant to the dictionary, associating it with the sample
            variants[variant_key][sample].add(1)  # presence is marked as 1
            samples.add(sample)

    return variants, sorted(samples)

def build_presence_absence_matrix(variants, samples, key_columns, output_file):
    """Build and print the presence/absence matrix for variants across samples."""
    # Open the output file in write mode
    with open(output_file, 'w') as f:
        # Write the header with sample names
        f.write("\t".join(key_columns) + "\t" + "\t".join(samples) + "\n")

        # Iterate through each variant and output the presence/absence for each sample
        for variant_key, sample_presence in variants.items():
            row = list(variant_key)

            # Fill in 1 or 0 for each sample
            row += [str(1 if sample in sample_presence else 0) for sample in samples]

            # Write the row as a tab-separated string
            f.write("\t".join(row) + "\n")

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Build presence/absence matrix based on key columns and sample column given by the user.")
    parser.add_argument('-i', '--input_file', required=True, help="Input TSV file path")
    parser.add_argument('-k', '--key_columns', nargs='+', required=True, help="List of key columns to use for merging")
    parser.add_argument('-s', '--sample_column', required=True, help="Name of the column indicated the sample")
    parser.add_argument('-o', '--output_file', required=True, help="Output TSV file path")

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    key_columns = args.key_columns
    sample_column = args.sample_column
    output_file = args.output_file

    # Parse the input file and extract variants based on the key columns
    variants, samples = parse_input_file(input_file, key_columns, sample_column)

    # Build and output the presence/absence matrix
    build_presence_absence_matrix(variants, samples, key_columns, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
