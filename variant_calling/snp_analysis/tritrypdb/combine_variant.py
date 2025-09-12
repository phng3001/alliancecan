# P=NP
import sys
import argparse

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
    except ImportError:
        print("Pandas is not available. Please install/load it.")
        sys.exit(1)

def combine_variant_rows(input_file, output_file):
    # Read the input TSV file
    df = pd.read_csv(input_file, sep='\t', dtype=str)

    # Define the columns to group by
    group_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'FORMAT', 'SAMPLE_INFO', 'SAMPLE']

    # Define the columns to merge (all columns except the group columns)
    merge_cols = [col for col in df.columns if col not in group_cols]

    # Group by the specified columns and combine the other columns with '|'
    def merge_values(group):
        # For each column that needs to be merged, join non-identical values with '|'
        merged_row = group.iloc[0].copy()
        for col in merge_cols:
            unique_values = group[col].dropna().unique()
            if len(unique_values) == 0:
                merged_row[col] = ''  # Set as empty string if no values are present
            elif len(unique_values) == 1:
                merged_row[col] = unique_values[0]
            else:
                merged_row[col] = '|'.join(unique_values)
        return merged_row

    # Apply the merging function to each group
    combined_df = df.groupby(group_cols).apply(merge_values).reset_index(drop=True)

    # Write the output to a new TSV file
    combined_df.to_csv(output_file, sep='\t', index=False)

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Combine variant rows based on specified columns. Correct the columns if needed inside the script.")
    parser.add_argument('-i', '--input', required=True, help="Input TSV file path")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file path")

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input
    output_file = args.output

    # Check if Pandas is available
    check_pandas()

    # Run the function with provided arguments
    combine_variant_rows(input_file, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
