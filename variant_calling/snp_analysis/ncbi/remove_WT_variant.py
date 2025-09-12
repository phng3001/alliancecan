# P=NP
import sys
import argparse

def parse_vcf(file_path):
    """Parse a VCF file and returns a set of variants as tuples (CHROM, POS, REF, ALT)."""
    variants = set()
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("#"):  # Skip header lines
                continue
            columns = line.strip().split('\t')
            chrom = columns[0]
            pos = int(columns[1])
            ref = columns[3]
            alt = columns[4]
            variants.add((chrom, pos, ref, alt))
    return variants

def contains_stretch_of_N(sequence):
    """Check if a sequence contains a stretch of 'N'."""
    return 'N' in sequence

def filter_unique_variants(wild_type_vcf, mutant_vcf, output_vcf):
    """Compare the wild type and mutant VCFs and outputs unique variants from the mutant.
    Also remove N-stretch-containing variants."""
    # Parse the VCF files
    wild_type_variants = parse_vcf(wild_type_vcf)

    # Open the output VCF file for writing
    with open(output_vcf, 'w') as out_file:
        # Write the VCF header
        with open(mutant_vcf, 'r') as mutant_file:
            for line in mutant_file:
                if line.startswith("#"):
                    out_file.write(line)  # Copy the header from mutant VCF

        # Read the mutant VCF and write only unique variants
        with open(mutant_vcf, 'r') as mutant_file:
            for line in mutant_file:
                if line.startswith("#"):  # Skip header lines
                    continue
                columns = line.strip().split('\t')
                chrom = columns[0]
                pos = int(columns[1])
                ref = columns[3]
                alt = columns[4]
                
                # Skip lines with a stretch of N in REF or ALT
                if contains_stretch_of_N(ref) or contains_stretch_of_N(alt):
                    continue
                
                # Check if the variant is unique to the mutant and write to the output VCF
                if (chrom, pos, ref, alt) not in wild_type_variants:
                    out_file.write(line)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description="Remove wild type variants and N-stretch-containing variants from mutant VCF file. Output mutant unique variants.")
    parser.add_argument('-w', '--wild_type', required=True, help="Wild type VCF file")
    parser.add_argument('-m', '--mutant', required=True, help="Mutant VCF file")
    parser.add_argument('-o', '--output', required=True, help="Output VCF file name")

    # Parse the arguments
    args = parser.parse_args()
    wild_type_vcf = args.wild_type
    mutant_vcf = args.mutant
    output_vcf = args.output

    # Call the filtering function with the provided arguments
    filter_unique_variants(wild_type_vcf, mutant_vcf, output_vcf)



if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()