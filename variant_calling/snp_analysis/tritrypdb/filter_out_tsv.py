# P=NP
import sys
import argparse

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
        print(f"Pandas is available. Version: {pd.__version__}", file=sys.stderr)
    except ImportError:
        print("Pandas is not available. Please install/load it.", file=sys.stderr)
        sys.exit(1)

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Remove rows from a TSV file where the specified column matches a given value.")

    parser.add_argument('-i', '--input', required=True, help="Input TSV file")
    parser.add_argument('-c', '--column', required=True, help="Column name to filter on")
    parser.add_argument('-f', '--filter', required=True, help="Value to filter out (i.e., exclude rows with this value)")
    parser.add_argument('-o', '--output', required=True, help="Output TSV file")

    args = parser.parse_args()

    # Check if Pandas is available
    check_pandas()

    # Load the TSV file
    try:
        df = pd.read_csv(args.input, sep='\t')
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)

    # Validate column
    if args.column not in df.columns:
        print(f"Error: Column '{args.column}' not found in the input file.")
        sys.exit(1)

    # Perform the filtering
    filtered_df = df[df[args.column] != args.filter]

    # Save to output file
    try:
        filtered_df.to_csv(args.output, sep='\t', index=False)
        print(f"Filtered data saved to '{args.output}'")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)



if __name__ == "__main__":
    main()
