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

def handle_featureCounts(input_file, output_file):
    # Load featureCounts input file (skip comment lines)
    df = pd.read_csv(input_file, sep='\t', comment='#')

    # Rename the column Geneid to transcript_id
    df = df.rename(columns={'Geneid': 'transcript_id'})

    # Set transcript_id column as index
    df.set_index('transcript_id', inplace=True)

    # Drop non-count columns (Chr, Start, End, Strand, Length)
    count_data = df.drop(columns=['Chr', 'Start', 'End', 'Strand', 'Length'])

    # Remove the ".bam" extension from column names
    count_data.columns = count_data.columns.str.replace('.bam', '', regex=False)

    # Save to output
    count_data.to_csv(output_file, sep=',')
    print(f"Handled csv file saved to: {output_file}")



def main():
    parser = argparse.ArgumentParser(description="Convert featureCounts tab delimited output to csv. Only id and count columns will be kept, id column will be named 'transcript_id'. The .bam extension in sample names will be removed.")
    parser.add_argument("-i", "--input", required=True, help="Input featureCounts file")
    parser.add_argument("-o", "--output", required=True, help="Output csv file")
    args = parser.parse_args()

    check_pandas()

    handle_featureCounts(args.input, args.output)



if __name__ == "__main__":
    main()
