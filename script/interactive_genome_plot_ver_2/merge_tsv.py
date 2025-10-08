# P=NP
import glob
import os
import argparse
import sys
import re

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
        print(f"Pandas is available. Version: {pd.__version__}", file=sys.stderr)
    except ImportError:
        print("Pandas is not available. Please install/load it.", file=sys.stderr)
        sys.exit(1)

def chrom_sort_key(chrom):
    """
    Tries to produce a natural sort key for chromosome names.
    For example: chr1, chr2, ..., chr10, chrX, chrY
    """
    if isinstance(chrom, str):
        # Extract number if possible
        match = re.match(r"^chr(\d+)$", chrom)
        if match:
            return (0, int(match.group(1)))
        # Put chrX, chrY, chrM after numbers
        match = re.match(r"^chr([A-Za-z]+)$", chrom)
        if match:
            order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
            return (1, order.get(match.group(1).upper(), 100))
        return (2, chrom)
    return (3, chrom)

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple TSV files on user-defined common columns "
            "and extract selected target columns per sample."
        )
    )
    parser.add_argument(
        "-i", "--inputs",
        required=True,
        help="Glob pattern (e.g., '*.tsv') or space-separated list of TSV input files"
    )
    parser.add_argument(
        "-c", "--common",
        required=True,
        help="Comma-separated list of columns common to all files (e.g., '#CHROM,START,END')"
    )
    parser.add_argument(
        "-t", "--target",
        required=True,
        help="Comma-separated list of columns to extract from each file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV filename"
    )
    parser.add_argument(
        "--sort",
        action="store_true",
        help="Sort output by the common columns."
    )

    args = parser.parse_args()

    check_pandas()

    # Expand file list
    if "*" in args.inputs or "?" in args.inputs:
        input_files = glob.glob(args.inputs)
    else:
        input_files = args.inputs.split()

    if not input_files:
        print("❌ No input files found.", file=sys.stderr)
        sys.exit(1)

    common_cols = [c.strip() for c in args.common.split(",")]
    target_cols = [c.strip() for c in args.target.split(",")]

    merged_df = None

    for filepath in input_files:
        sample_name = os.path.splitext(os.path.basename(filepath))[0]
        print(f"Processing {sample_name}...")

        df = pd.read_csv(filepath, sep="\t", dtype=str)

        # Check columns
        missing_common = [c for c in common_cols if c not in df.columns]
        missing_target = [c for c in target_cols if c not in df.columns]
        if missing_common:
            print(f"⚠️ Skipping {sample_name}: missing common columns {missing_common}")
            continue
        if missing_target:
            print(f"⚠️ Skipping {sample_name}: missing target columns {missing_target}")
            continue

        # Keep only common + target columns
        keep_cols = common_cols + target_cols
        df = df[keep_cols].copy()

        # Rename target columns with sample prefix if multiple target selected
        rename_map = {
            col: f"{sample_name}_{col}" if len(target_cols) > 1 else sample_name
            for col in target_cols
        }
        df.rename(columns=rename_map, inplace=True)

        # Merge
        if merged_df is None:
            merged_df = df
        else:
            merged_df = pd.merge(merged_df, df, on=common_cols, how="outer")

    if merged_df is None:
        print("❌ No valid data merged.")
        sys.exit(1)

    # Optional sorting
    if args.sort:
        print("Sorting by common columns...")
        if "#CHROM" in common_cols:
            # Handle chromosome natural sort
            merged_df["#CHROM_sort"] = merged_df["#CHROM"].apply(chrom_sort_key)
            sort_cols = ["#CHROM_sort"] + [c for c in common_cols if c != "#CHROM"]
            merged_df[sort_cols[1:]] = merged_df[sort_cols[1:]].apply(pd.to_numeric, errors="ignore")
            merged_df = merged_df.sort_values(by=sort_cols, ascending=True).drop(columns=["#CHROM_sort"])
        else:
            merged_df[common_cols] = merged_df[common_cols].apply(pd.to_numeric, errors="ignore")
            merged_df = merged_df.sort_values(by=common_cols, ascending=True)

    merged_df.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Merged file written to: {args.output}")



if __name__ == "__main__":
    main()