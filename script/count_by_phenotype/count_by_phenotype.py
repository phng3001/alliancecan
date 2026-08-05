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

def process_data(presence_absence_file, phenotype_file, output_file):
    # Load presence-absence matrix
    presence_absence_df = pd.read_csv(presence_absence_file, sep='\t')

    # Set the first column (variant IDs) as the index for presence_absence_df
    presence_absence_df.set_index(presence_absence_df.columns[0], inplace=True)

    # Transpose presence-absence matrix to have strains as index
    presence_absence_df = presence_absence_df.T

    # Load phenotype data
    phenotype_df = pd.read_csv(phenotype_file, sep='\t')

    # Set strain as the index for phenotype_df for merging
    phenotype_df.set_index(phenotype_df.columns[0], inplace=True)

    # Merge presence-absence matrix with phenotype table
    merged_df = presence_absence_df.merge(phenotype_df, left_index=True, right_index=True)

    # Extract phenotype column name
    phenotype_col = phenotype_df.columns[0]

    # Group by phenotype and sum the presence-absence values for each variant, then transpose to have variants as rows and phenotypes as columns
    var_counts_df = merged_df.groupby(phenotype_col).sum().T

    # Save the result
    var_counts_df.to_csv(output_file, header=True, index=True, index_label=phenotype_col, sep='\t')
    print(f"Results saved to {output_file}")

def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Process variant presence-absence matrix and phenotype data to count variant per phenotype.
        """
        )
    
    parser.add_argument('--pres', required=True, help='Presence-absence matrix TSV file')
    parser.add_argument('--phenotypes', required=True, help='Phenotype data TSV file')
    parser.add_argument('--output', required=True, help='Output TSV file')

    # Parse arguments
    args = parser.parse_args()

    check_pandas()

    # Process data
    process_data(args.pres, args.phenotypes, args.output)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
