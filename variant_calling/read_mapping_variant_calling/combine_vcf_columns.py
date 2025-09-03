# P=NP
import sys

def combine_columns(vcf_line):
    columns = vcf_line.strip().split("\t")

    # Extract the fields related to genotype from both columns
    snp_data = columns[9]  # SNPs column
    indel_data = columns[10]  # INDELs column

    # Check if either column has missing data and handle accordingly
    if snp_data == "./.:.:.:.:.":
        combined_data = indel_data
    elif indel_data == "./.:.:.:.:.":
        combined_data = snp_data
    else:
        # Combine both fields - you can customize how to handle combined data
        combined_data = f"{snp_data}/{indel_data}"  # Simple concatenation, adjust as needed

    # Replace the two columns with the combined one
    columns[9] = combined_data
    del columns[10]  # Remove the indel column

    # Reconstruct the line
    return "\t".join(columns)

def process_vcf(input_vcf, output_vcf, new_header):
    # Open the input VCF and apply the combination to each line
    with open(input_vcf, "r") as infile, open(output_vcf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                # Modify the header line
                if line.startswith("#CHROM"):
                    # Create new header line
                    header = f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{new_header}\n"
                    outfile.write(header)
                else:
                    outfile.write(line)  # Write other header lines unchanged
            else:
                outfile.write(combine_columns(line) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python combine_vcf_columns.py <input.vcf> <output.vcf> <new_header>")
        sys.exit(1)

    input_vcf = sys.argv[1]
    output_vcf = sys.argv[2]
    new_header = sys.argv[3]  # sample name

    process_vcf(input_vcf, output_vcf, new_header)

