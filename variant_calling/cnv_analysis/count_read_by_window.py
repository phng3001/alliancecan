# P=NP
import sys
import subprocess
import argparse
import os

def check_samtools():
    """Check if Samtools is available."""
    try:
        subprocess.run(["samtools", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    except FileNotFoundError:
        print("Samtools is not available. Please install/load it.")
        exit(1)

def index_files(bam_file, fasta_file):
    """Index BAM and FASTA files if not already indexed."""
    # Check and index the BAM file
    bam_index = f"{bam_file}.bai"
    if not os.path.exists(bam_index):
        print(f"Indexing BAM file: {bam_file}")
        subprocess.run(["samtools", "index", bam_file], check=True)

    # Check and index the FASTA file
    fasta_index = f"{fasta_file}.fai"
    if not os.path.exists(fasta_index):
        print(f"Indexing FASTA file: {fasta_file}")
        subprocess.run(["samtools", "faidx", fasta_file], check=True)

def count_reads_by_window(bam_file, fasta_file, window_size, output_file):
    """Count reads from BAM file for genomic windows of given size."""

    # Ensure necessary files are indexed
    index_files(bam_file, fasta_file)

    # Parse genome fasta index (fai) file to get reference lengths
    genome_fai = f"{fasta_file}.fai"
    genome_windows = []
    with open(genome_fai, "r") as fai:
        for line in fai:
            chrom, length = line.split()[:2]
            length = int(length)
            for start in range(0, length, window_size):
                end = min(start + window_size, length)
                genome_windows.append((chrom, start, end))

    # Prepare output file
    with open(output_file, "w") as out:
        out.write("#CHROM\tSTART\tEND\tNB_READS\n")

        for chrom, start, end in genome_windows:
            if start == 0:
                print(f"Start processing chromosome {chrom} in {bam_file}")
            # print(f"Start processing {chrom} {start}-{end} region")

            # Use samtools to count reads in the current window
            region = f"{chrom}:{start}-{end}"
            cmd = ["samtools", "view", "-c", bam_file, region]
            try:
                result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
                nb_reads = result.stdout.strip()
                out.write(f"{chrom}\t{start}\t{end}\t{nb_reads}\n")
            except subprocess.CalledProcessError as e:
                print(f"Error counting reads in {region}: {e.stderr}")
                out.write(f"{chrom}\t{start}\t{end}\t0\n")  # Write 0 read for problematic regions

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Count reads in BAM file by genomic windows.")
    parser.add_argument('-i', '--bam_file', required=True, help="Input BAM file path")
    parser.add_argument('-r', '--fasta_file', required=True, help="Reference FASTA file path")
    parser.add_argument('-w', '--window_size', required=True, type=int, help="Size of each genomic window (e.g., 5000).")
    parser.add_argument('-o', '--output_file', required=True, help="Output TSV file path")

    # Parse the arguments
    args = parser.parse_args()
    bam_file = args.bam_file
    fasta_file = args.fasta_file
    window_size = args.window_size
    output_file = args.output_file

    # Check if Samtools is available
    check_samtools()

    # Run the function with provided arguments
    count_reads_by_window(bam_file, fasta_file, window_size, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
    