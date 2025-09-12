# P=NP
import argparse
import os
import sys

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
    except ImportError:
        print("Pandas is not available. Please install/load it.")
        sys.exit(1)

def parse_gff(gff_file, feature_type):
    """Parse the GFF file and extract information based on the given feature type."""
    genes = []
    missing_description_count = 0  # Count locus tags missing product/description
    total_genes = 0  # Track total number of relevant genes

    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            chrom, source, feature, start, end, score, strand, frame, attribute = fields
            
            if feature == feature_type:
                total_genes += 1
                attribute_dict = {}
                for attr in attribute.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attribute_dict[key] = value

                locus_tag = attribute_dict.get('ID') or attribute_dict.get('locus_tag', '-')
                gene_name = attribute_dict.get('Name') or attribute_dict.get('gene', '-')
                description = attribute_dict.get('product') or attribute_dict.get('description', None)

                if description is None:
                    missing_description_count += 1
                    description = '-'  # Default to '-' if not found

                genes.append({
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'locus_tag': locus_tag,
                    'gene_name': gene_name,
                    'description': description
                })

    # Print a warning only if all locus tags are missing product/description
    if total_genes > 0 and missing_description_count == total_genes:
        print(f"Warning: No 'product' or 'description' found for any of the {total_genes} locus tags in the GFF file.")

    return genes

def find_genes_in_region(chrom, start, end, genes):
    """Find all genes within a given region."""
    region_genes = []
    for gene in genes:
        if gene['chrom'] == chrom and gene['start'] <= end and gene['end'] >= start:
            region_genes.append(
                f"({gene['locus_tag']},{gene['gene_name']},{gene['description']})"
            )
    return ';'.join(region_genes)

def annotate_genomic_region(input_file, gff_file, feature_type, output_file):
    """Annotate the input TSV file with gene information from the GFF file."""
    # Load the input TSV file
    df = pd.read_csv(input_file, sep='\t')

    # Ensure required columns are present
    if not {'#CHROM', 'START', 'END'}.issubset(df.columns):
        print("Error: Input TSV file must contain '#CHROM', 'START', and 'END' columns.")
        sys.exit(1)

    # Parse the GFF file
    genes = parse_gff(gff_file, feature_type)

    # Annotate each row in the TSV file
    # df['GENE'] = df.apply(
    #     lambda row: find_genes_in_region(row['#CHROM'], row['START'], row['END'], genes), 
    #     )
    
    df['GENE'] = df.apply(
        lambda row: find_genes_in_region(row['#CHROM'], row['START'], row['END'], genes) 
        if find_genes_in_region(row['#CHROM'], row['START'], row['END'], genes) else '-', axis=1
        )

    # Save the annotated TSV file
    df.to_csv(output_file, sep='\t', index=False)

def main():
    # Define command line arguments
    parser = argparse.ArgumentParser(description="""Annotate genomic regions with gene information retrieved from a GFF file.
                                     Input file must contain these columns: '#CHROM', 'START', and 'END'.
                                     Gene information will be written to a new column named 'GENE'.""")
    parser.add_argument('-i', '--input_file', required=True, help="Input TSV file path")
    parser.add_argument('-r', '--gff_file', required=True, help="Input GFF file path")
    parser.add_argument('-f', '--feature', default='gene', help="Feature type to annotate (e.g., 'gene', 'CDS'). Default is 'gene'")
    parser.add_argument('-o', '--output_file', required=True, help="Output TSV file path")
    
    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    gff_file = args.gff_file
    feature_type = args.feature
    output_file = args.output_file
    
    # Check if Pandas is available
    check_pandas()

    # Annotate the TSV file
    annotate_genomic_region(input_file, gff_file, feature_type, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
