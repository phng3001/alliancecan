# P=NP
import sys
import argparse
import csv

def process_data(args):

    # Read second file into a dictionary
    lookup = {}

    with open(args.file2, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        if args.key2 not in reader.fieldnames:
            raise ValueError(f"Column '{args.key2}' not found in {args.file2}")

        file2_columns = [c for c in reader.fieldnames if c != args.key2]

        for row in reader:
            lookup[row[args.key2]] = row

    # Merge
    with open(args.file1, newline="") as fin, \
         open(args.output, "w", newline="") as fout:

        reader = csv.DictReader(fin, delimiter="\t")

        if args.key1 not in reader.fieldnames:
            raise ValueError(f"Column '{args.key1}' not found in {args.file1}")

        output_columns = reader.fieldnames + file2_columns

        writer = csv.DictWriter(
            fout,
            fieldnames=output_columns,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in reader:
            merged = row.copy()

            match = lookup.get(row[args.key1])

            if match:
                for col in file2_columns:
                    merged[col] = match[col]
            else:
                for col in file2_columns:
                    merged[col] = ""

            writer.writerow(merged)

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Merge two TSV files using two key columns.
        The first file is used as the base (left join).
        """
    )

    parser.add_argument("--file1", required=True, help="First TSV file (rows to keep)")
    parser.add_argument("--file2", required=True, help="Second TSV file")
    parser.add_argument("--key1", required=True, help="Key column in the first file")
    parser.add_argument("--key2", required=True, help="Key column in the second file")
    parser.add_argument("--output", required=True, help="Output TSV file")

    # Parse arguments
    args = parser.parse_args()

    # Process data
    process_data(args)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
