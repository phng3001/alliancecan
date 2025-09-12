# P=NP
import sys

def parse_vcf(vcf_file, sample_name):
    """Parse VCF file and extract variant information."""
    variants = []
    
    with open(vcf_file, 'r') as file:
        for line in file:
            if line.startswith("#"):
                continue   # Skip header lines

            # Handle VCF's fields           
            columns = line.strip().split('\t')
            chrom = columns[0]        # CHROM
            pos = columns[1]          # POS
            ref = columns[3]          # REF
            alt = columns[4]          # ALT
            qual = columns[5]         # QUAL
            filter = columns[6]       # FILTER
            info = columns[7]         # INFO
            format_field = columns[8] # FORMAT
            sample_info = columns[9]  # Sample information

            # Extracting GENOTYPE (GT)
            format_keys = format_field.split(':')
            sample_keys = sample_info.split(':')
            genotype = sample_keys[format_keys.index('GT')] if 'GT' in format_keys else "NA"

            # Find the ANN field in the INFO column
            ann_field = [field for field in info.split(';') if field.startswith("ANN=")]
            if not ann_field:
                continue

            annotations = ann_field[0][4:].split(',')
            # Handle ANN's fields
            for annotation in annotations:
                fields = annotation.split('|')
                allele = fields[0]
                annotation_type = fields[1]
                annotation_impact = fields[2]
                gene_name = fields[3]
                gene_id = fields[4]
                feature_type = fields[5]
                feature_id = fields[6]
                transcript_bioType = fields[7]
                rank = fields[8]
                nu_change = fields[9]
                aa_change = fields[10]
                cdna_pos = fields[11]
                cds_pos = fields[12]
                aa_pos = fields[13]
                distance = fields[14]
                notes = fields[15]

                # Store the variant information
                variant_info = (
                    chrom, pos, ref, alt, qual, format_field, sample_info, genotype, gene_id,
                    feature_id, annotation_type, nu_change, aa_change, cds_pos, aa_pos, sample_name
                )
                variants.append(variant_info)

    return variants

def aggregate_variants(sample_vcf_files):
    """Output unique variants from the given VCF files."""
    unique_variants = set()   # Use a set to store unique variants

    # Process each VCF file and collect unique variants
    for sample in sample_vcf_files:
        sample_name = sample.split('/')[-1]   # Extract sample name from the file path
        variants = parse_vcf(sample, sample_name)
        for variant in variants:
            unique_variants.add(variant)   # Add to set to ensure uniqueness

    # Print the header
    print("#CHROM\tPOS\tREF\tALT\tQUAL\tFORMAT\tSAMPLE_INFO\tGENOTYPE\tGENE_ID\tFEATURE_ID\tVARIANT_TYPE\tNU_CHANGE\tAA_CHANGE\tCDS_POS/CDS_LENGTH\tAA_POS/AA_LENGTH\tSAMPLE")

    # Print all unique variants
    for variant in sorted(unique_variants):
        print('\t'.join(variant))



if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python aggregate_variant.py <sample1.vcf> <sample2.vcf> ...")
        sys.exit(1)
    
    sample_vcfs = sys.argv[1:]
    aggregate_variants(sample_vcfs)