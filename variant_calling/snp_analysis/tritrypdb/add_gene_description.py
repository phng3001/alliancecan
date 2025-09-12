# P=NP
import sys
import argparse
import re

def check_pandas():
    """Check if Pandas is available."""
    try:
        global pd
        import pandas as pd
    except ImportError:
        print("Pandas is not available. Please install/load it.")
        sys.exit(1)

def parse_gff(gff_file):
    """Parse GFF file to extract annotations and gene names."""
    annotations = {}
    gene_names = {}

    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith("#"):
                fields = line.split('\t')
                attributes = fields[8]
                gene_id = None
                gene_name = None
                description = None

                # TriTrypDB gff
                for attribute in attributes.split(';'):
                    if attribute.startswith('ID='):
                        gene_id = attribute.split('=')[1].strip()
                    elif attribute.startswith('Name='):
                        gene_name = attribute.split('=')[1].strip()
                    elif attribute.startswith('description='):
                        description = attribute.split('=')[1].strip()

                if gene_id:
                    if description:
                        annotations[gene_id] = description
                    if gene_name:
                        gene_names[gene_id] = gene_name

    return annotations, gene_names

def extract_gene_ids(gene_id):
    """Clean and extract gene IDs."""
    # Remove "GENE_", "exon_", "-p*", "-CDS*", "-E*", "_gene", "gene-" if existed
    gene_ids = re.sub(r'GENE_|exon_|-p\d*-CDS\d*|-E\d*|_gene|gene-', '', gene_id)
    return gene_ids

def add_gene_description(input_file, gff_annotations, gff_gene_names, output_file):
    """Add gene names and descriptions to the input file."""
    # Load input file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t', dtype=str)
    
    # Check if 'GENE_ID' column exists
    if 'GENE_ID' not in df.columns:
        raise ValueError("The input file is missing the required column: 'GENE_ID'")
    
    # Add 'GENE_NAME' and 'DESCRIPTION' columns
    df['GENE_NAME'] = "-"
    df['DESCRIPTION'] = "-"
    
    # Process each row to add annotations
    for idx, row in df.iterrows():
        gene_id = row['GENE_ID']
        gene_ids = extract_gene_ids(gene_id)

        if "-" in gene_ids: # Intergenic
            gene_id_list = gene_ids.split("-")
            descriptions = [gff_annotations.get(gid, "-") for gid in gene_id_list]
            gene_names = [gff_gene_names.get(gid, "-") for gid in gene_id_list]
            df.at[idx, 'DESCRIPTION'] = "|".join(descriptions)
            df.at[idx, 'GENE_NAME'] = "|".join(gene_names)
        else:
            df.at[idx, 'DESCRIPTION'] = gff_annotations.get(gene_ids, "-")
            df.at[idx, 'GENE_NAME'] = gff_gene_names.get(gene_ids, "-")
    
    # Reorder columns to place 'GENE_NAME' and 'DESCRIPTION' after 'GENE_ID'
    gene_id_index = df.columns.get_loc('GENE_ID')
    columns = list(df.columns)
    columns.remove('GENE_NAME')
    columns.remove('DESCRIPTION')
    columns = (
        columns[:gene_id_index + 1] +  # Columns before and including 'GENE_ID'
        ['GENE_NAME', 'DESCRIPTION'] +  # Add the new columns
        columns[gene_id_index + 1:]  # Columns after 'GENE_ID'
    )
    df = df[columns]
    
    # Save the updated DataFrame to the output file
    df.to_csv(output_file, sep='\t', index=False)

def main():
    # Define command-line arguments
    parser = argparse.ArgumentParser(description="Add gene descriptions and names to input file based on TriTrypDB GFF file.")
    parser.add_argument('-i', '--input_file', required=True, help="Input file")
    parser.add_argument('-g', '--gff_file', required=True, help="GFF file")
    parser.add_argument('-o', '--output_file', required=True, help="Output file name")

    # Parse the arguments
    args = parser.parse_args()
    input_file = args.input_file
    gff_file = args.gff_file
    output_file = args.output_file

    # Check if Pandas is available
    check_pandas()

    # Parse the GFF file
    gff_annotations, gff_gene_names = parse_gff(gff_file)

    # Add gene descriptions and names
    add_gene_description(input_file, gff_annotations, gff_gene_names, output_file)

if __name__ == "__main__":
    # Check if no arguments are provided
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
