# P=NP
import argparse
import sys
import pandas as pd

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
        print(f"Pandas is available. Version: {pd.__version__}", file=sys.stderr)
    except ImportError:
        print("Pandas is not available. Please install/load it.", file=sys.stderr)
        sys.exit(1)

def process_data(pres_file, sample_file, start_col, output_file):
    """Subset a presence/absence matrix by samples and remove all-zero rows."""

    # Read presence/absence matrix
    try:
        df = pd.read_csv(pres_file, sep="\t", dtype=str)
    except Exception as e:
        sys.exit(f"Error reading presence/absence matrix: {e}")

    # Check starting column
    if start_col < 1:
        sys.exit("Error: start_col must be >= 1.")

    start_idx = start_col - 1

    if start_idx >= len(df.columns):
        sys.exit(
            f"Error: start_col ({start_col}) exceeds the number of columns in the matrix ({len(df.columns)})."
        )

    # Read sample list
    try:
        with open(sample_file, "r") as f:
            samples = [line.strip() for line in f if line.strip()]
    except Exception as e:
        sys.exit(f"Error reading sample list: {e}")

    if not samples:
        sys.exit("Error: Sample list is empty.")

    # Check for duplicate sample names
    if len(samples) != len(set(samples)):
        sys.exit("Error: Duplicate sample names found in sample list.")

    # Identify sample columns
    sample_columns = list(df.columns[start_idx:])
    description_columns = list(df.columns[:start_idx])

    # Check whether requested samples exist
    missing = [sample for sample in samples if sample not in sample_columns]

    if missing:
        sys.exit(
            "Error: The following samples were not found in the matrix:\n"
            + "\n".join(missing)
        )

    # Subset columns (preserve sample list order)
    selected_columns = description_columns + samples
    df_subset = df[selected_columns].copy()

    # Check and convert presence/absence values
    try:
        values = df_subset[samples].apply(pd.to_numeric, errors="raise")
    except Exception as e:
        sys.exit(f"Error: Non-numeric values found in sample columns: {e}")

    # Remove rows where all selected samples have value 0
    df_subset = df_subset.loc[values.ne(0).any(axis=1)].copy()

    # Write output
    try:
        df_subset.to_csv(output_file, sep="\t", index=False)
    except Exception as e:
        sys.exit(f"Error writing output file: {e}")

    # Summary
    print(f"Input matrix: {pres_file}")
    print(f"Total samples in matrix: {len(sample_columns)}")
    print(f"Selected samples: {len(samples)}")
    print(f"Input features: {len(df)}")
    print(f"Remaining features: {len(df_subset)}")
    print(f"Removed features (all zero): {len(df) - len(df_subset)}")
    print(f"Output written to: {output_file}")



def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        description="Subset a presence/absence matrix by samples (columns)."
    )

    parser.add_argument("--pres", required=True, help="Presence/absence matrix (TSV)")
    parser.add_argument("--sample_list", required=True, help="Sample list (one per line)")
    parser.add_argument("--start_col", type=int, default=2, help="1-based column number where sample columns start (default: 2)")
    parser.add_argument("--output", required=True, help="Output subset matrix (TSV)")

    # Parse arguments
    args = parser.parse_args()

    # Check pandas
    check_pandas()

    # Process data
    process_data(args.pres, args.sample_list, args.start_col, args.output)



if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
