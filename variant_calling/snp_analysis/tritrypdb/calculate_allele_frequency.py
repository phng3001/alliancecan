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

def calculate_allele_freq(row):
    # Split SAMPLE_INFO to get AD and DP
    sample_info = row['SAMPLE_INFO']
    format_info = row['FORMAT']

    # Check if AD and DP are present in the FORMAT column
    format_list = format_info.split(':')    
    try:
        ad_index = format_list.index('AD')  # Find the index of AD
        dp_index = format_list.index('DP')  # Find the index of DP
    except ValueError:
        # If AD or DP is missing from FORMAT, return an error message
        return "ERROR: AD or DP not found in FORMAT"

    # Extract the values for AD and DP from SAMPLE_INFO
    sample_info_list = sample_info.split(':')
    try:
        ad = sample_info_list[ad_index]
        dp = sample_info_list[dp_index]
    except IndexError:
        # If SAMPLE_INFO is improperly formatted
        return "ERROR: Invalid SAMPLE_INFO format"

    # Convert AD and DP to integers
    try:
        ad_values = list(map(int, ad.split(',')))  # AD values are separated by commas
        dp_value = int(dp)  # DP is a single integer
    except ValueError:
        return "ERROR: Invalid AD or DP values"

    # Calculate allele frequencies for each AD value
    if dp_value > 0:
        allele_freqs = [(ad_val / dp_value) for ad_val in ad_values]
    else:
        allele_freqs = [0] * len(ad_values)  # If DP is 0, set allele frequencies to 0
    
    # Join the allele frequencies with commas
    return ','.join(f"{freq:.2f}" for freq in allele_freqs)

def calculate_allele_frequencies(input_file, output_file):
    # Read the input TSV file
    try:
        df = pd.read_csv(input_file, sep='\t', dtype=str)
    except Exception as e:
        print(f"Error: Unable to read input file '{input_file}'. {e}")
        sys.exit(1)

    # Check if required columns are present in the input file
    required_columns = ['SAMPLE_INFO', 'FORMAT']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: The following required columns are missing: {', '.join(missing_columns)}")
        sys.exit(1)

    # Apply the allele frequency calculation
    df['ALLELE_FREQ'] = df.apply(calculate_allele_freq, axis=1)

    # Reorder the columns to insert 'ALLELE_FREQ' after 'SAMPLE_INFO'
    cols = list(df.columns)
    # Find the position of 'SAMPLE_INFO' and insert 'ALLELE_FREQ' after it
    sample_info_index = cols.index('SAMPLE_INFO')
    cols.insert(sample_info_index + 1, cols.pop(cols.index('ALLELE_FREQ')))

    # Rearrange the columns
    df = df[cols]

    # Write the output to a new TSV file
    try:
        df.to_csv(output_file, sep='\t', index=False)
    except Exception as e:
        print(f"Error: Unable to write output file '{output_file}'. {e}")
        sys.exit(1)

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Calculate allele frequencies from FORMAT and SAMPLE_INFO columns.")
    parser.add_argument('-i', '--input_file', required=True, help="Input TSV file path")
    parser.add_argument('-o', '--output_file', required=True, help="Output TSV file path")

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file

    # Check if Pandas is available
    check_pandas()

    # Run the function with provided arguments
    calculate_allele_frequencies(input_file, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
