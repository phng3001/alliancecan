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

def add_column_with_priority(input_file, ref_files, key_columns, column_name, column_position, output_file, verbose=False):
    # Load the input file into a DataFrame
    try:
        input_df = pd.read_csv(input_file, sep='\t', dtype=str)
    except Exception as e:
        print(f"Error: Unable to read input file '{input_file}'. {e}")
        sys.exit(1)

    input_df = input_df.copy()  # Avoid mutating the original input_df
    input_df[column_name] = None  # Initialize the user-specified column

    # Validate required key columns in the input file
    for col in key_columns:
        if col not in input_df.columns:
            raise ValueError(f"Input file is missing required column: {col}")

    # Process each reference file in order
    for ref_file in ref_files:
        try:
            ref_df = pd.read_csv(ref_file, sep='\t', dtype=str)
        except Exception as e:
            print(f"Error: Unable to read reference file '{ref_file}'. {e}")
            sys.exit(1)

        # Validate required key columns in the reference file
        for col in key_columns:
            if col not in ref_df.columns:
                raise ValueError(f"Reference file '{ref_file}' is missing required column: {col}")

        if column_name not in ref_df.columns:
            raise ValueError(f"The file '{ref_file}' must contain the '{column_name}' column.")

        # Check for duplicates in the reference file
        duplicates = ref_df.duplicated(subset=key_columns, keep=False)
        if duplicates.any():
            duplicate_rows = ref_df[duplicates].sort_values(by=key_columns)
            print(f"Warning: Duplicate rows found in reference file: {ref_file}")
            for idx, row in duplicate_rows.iterrows():
                print(f"Line {idx + 2}: {row.to_dict()}")

        # Merge input with the current ref file
        input_df = pd.merge(
            input_df,
            ref_df[key_columns + [column_name]],
            on=key_columns,
            how='left',
            suffixes=('', '_new')
        )

        # Update the user-specified column only where it's still missing
        input_df[column_name] = input_df[column_name].combine_first(input_df[f"{column_name}_new"])

        # Drop the temporary column
        input_df.drop(columns=[f"{column_name}_new"], inplace=True)

        if verbose:
            print(f"Reference file {ref_file} processed")
    
    # Move the column to the specified position if provided
    if column_position is not None:
        try:
            column_position = int(column_position)
            if column_position < 0 or column_position > len(input_df.columns):
                raise ValueError(f"Invalid column position: {column_position}")
            # Move the column to the specified position
            column_data = input_df.pop(column_name)
            input_df.insert(column_position, column_name, column_data)
        except Exception as e:
            print(f"Error: Unable to set column position. {e}")
            sys.exit(1)

    # Save the output to a new TSV file
    try:
        input_df.to_csv(output_file, sep='\t', index=False)
    except Exception as e:
        print(f"Error: Unable to write output file '{output_file}'. {e}")
        sys.exit(1)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description="Add a specified column from reference files to input file with priority order.")
    parser.add_argument('-i', '--input_file', required=True, help="Input TSV file path")
    parser.add_argument('-r', '--ref_files', nargs='+', required=True, help="List of reference files in priority order")
    parser.add_argument('-k', '--key_columns', nargs='+', required=True, help="List of key columns to use for merging")
    parser.add_argument('-c', '--column_name', required=True, help="Name of the column to add")
    parser.add_argument('-p', '--column_position', type=int, help="Zero-based position to place the new column")
    parser.add_argument('-o', '--output_file', required=True, help="Output TSV file path")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose mode")

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    ref_files = args.ref_files
    key_columns = args.key_columns
    column_name = args.column_name
    column_position = args.column_position
    output_file = args.output_file
    verbose = args.verbose

    # Check if Pandas is available
    check_pandas()

    # Run the function with provided arguments
    add_column_with_priority(input_file, ref_files, key_columns, column_name, column_position, output_file, verbose)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
      