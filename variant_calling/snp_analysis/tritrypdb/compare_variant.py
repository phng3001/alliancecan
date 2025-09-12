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

def load_variants(file_path):
    """Load the variants from the given file into a DataFrame."""
    try:
        df = pd.read_csv(file_path, sep='\t', header=0)
        df.columns = df.columns.str.strip()  # Strip whitespace from column names
        return df
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def compare_variants(dataset_dict, key_columns):
    """Compare multiple DataFrames and return a combined DataFrame with a CALLED_BY column."""
    # Concatenate all DataFrames with an additional SOURCE column
    combined_df = pd.concat(
        [df.assign(SOURCE=name)[key_columns + ['SOURCE']] for name, df in dataset_dict.items()],
        ignore_index=True
    ).fillna("-")  # Replace missing values with a placeholder

    # Aggregate sources for each unique combination of key columns
    result_df = (
        combined_df
        .groupby(key_columns, as_index=False)
        .agg(SOURCE=('SOURCE', lambda x: '|'.join(sorted(set(x)))))
    )
    
    # Rename SOURCE column to CALLED_BY
    result_df.rename(columns={'SOURCE': 'CALLED_BY'}, inplace=True)

    return result_df

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Compare variant files based on key columns.")
    parser.add_argument('-k', '--key_columns', nargs='+', required=True, help="List of key columns (space-separated)")
    parser.add_argument('-n', '--name', action='append', required=True, help="Algorithm name(s)")
    parser.add_argument('-f', '--file', action='append', required=True, help="Algorithm TSV data file(s)")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file name")

    args = parser.parse_args()

    # Check if Pandas is available
    check_pandas()

    # Validate input lengths
    if len(args.name) != len(args.file):
        print("Error: The number of algorithm names must match the number of data files.")
        sys.exit(1)

    # Load datasets into a dictionary
    dataset_dict = {
        name: load_variants(file_path) for name, file_path in zip(args.name, args.file)
    }

    # Perform comparison
    result_df = compare_variants(dataset_dict, args.key_columns)

    # Write results to output file
    try:
        result_df.to_csv(args.output, sep='\t', index=False)
    except Exception as e:
        print(f"Error writing to file {args.output}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Provide help message if no arguments are supplied
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
