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

def check_numpy():
    """Check if Numpy is available."""
    try:
        global np
        import numpy as np
    except ImportError:
        print("Numpy is not available. Please install/load it.")
        sys.exit(1)

def compute_log2_ratio(mutant_file, wt_file, output_file):
    # Load the mutant and wild type data
    mutant_df = pd.read_csv(mutant_file, sep="\t")
    wt_df = pd.read_csv(wt_file, sep="\t")

    # Verify the files have matching rows
    if not (mutant_df[["#CHROM", "START", "END"]].equals(wt_df[["#CHROM", "START", "END"]])):
        raise ValueError("The mutant and wild type files do not have matching genomic windows.")

    # Combine the data
    combined_df = mutant_df[["#CHROM", "START", "END"]].copy()
    combined_df["NB_READS_MUTANT"] = mutant_df["NB_READS"]
    combined_df["NB_READS_WT"] = wt_df["NB_READS"]

    # Calculate mutant/wt read count ratio for each window
    combined_df["NB_READS_RATIO_MUTANT/WT"] = combined_df["NB_READS_MUTANT"] / combined_df["NB_READS_WT"].replace(0, np.nan)

    # Calculate total read count
    total_reads_mutant = combined_df["NB_READS_MUTANT"].sum()
    total_reads_wt = combined_df["NB_READS_WT"].sum()

    # Calculate total mutant/wt read count ratio
    if total_reads_wt == 0:
        print("Warning: TOTAL_READS_WT is zero. Ratios will be set to NaN.")
        combined_df["TOTAL_READS_MUTANT"] = total_reads_mutant
        combined_df["TOTAL_READS_WT"] = np.nan
        combined_df["TOTAL_READS_RATIO_MUTANT/WT"] = np.nan
    else:
        combined_df["TOTAL_READS_MUTANT"] = total_reads_mutant
        combined_df["TOTAL_READS_WT"] = total_reads_wt
        combined_df["TOTAL_READS_RATIO_MUTANT/WT"] = total_reads_mutant / total_reads_wt

    # Calculate the log2 ratio of mutant/wt read counts for each window normalized by the total read count per sample
    combined_df["NORMALIZED_LOG2_READS_RATIO_MUTANT/WT"] = np.log2(
        combined_df["NB_READS_RATIO_MUTANT/WT"] / combined_df["TOTAL_READS_RATIO_MUTANT/WT"]
    ).replace([np.inf, -np.inf], np.nan)

    # Write the output
    combined_df.to_csv(output_file, sep="\t", index=False)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description="Compute the mutant/wild type read ratio in log2 for each corresponding genomic window. Read counts are normalized by the total read count per sample.")
    parser.add_argument('-i', '--mutant_file', required=True, help="Mutant read counts file path")
    parser.add_argument('-r', '--wt_file', required=True, help="Wild type read counts file path")
    parser.add_argument('-o', '--output_file', required=True, help="Output file path")
    
    # Parse the arguments
    args = parser.parse_args()
    mutant_file = args.mutant_file
    wt_file = args.wt_file
    output_file = args.output_file

    # Check if Pandas and Numpy are available
    check_pandas()
    check_numpy()

    # Run the function with provided arguments
    compute_log2_ratio(mutant_file, wt_file, output_file)



if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()