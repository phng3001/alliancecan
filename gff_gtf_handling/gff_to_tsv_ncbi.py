import sys
import argparse
import csv

def parse_gff3_attributes(attr_string):
    """Parse GFF attribute string (9th column of a GFF file) into a dictionary."""
    attr_dict = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attr_dict[key] = value
    return attr_dict

# NCBI Genbank GFF3 files
## Example 1:
### CDS
# CP103914.1	Genbank	gene	3908	4906	.	-	.	ID=gene-NXY56_000001;Name=NXY56_000001;gbkey=Gene;gene_biotype=protein_coding;locus_tag=NXY56_000001
# CP103914.1	Genbank	mRNA	3908	4906	.	-	.	ID=rna-LgyM4147_000005000.1;Parent=gene-NXY56_000001;gbkey=mRNA;locus_tag=NXY56_000001;orig_protein_id=gnl|NXY56|LgyM4147_000005000.1:CDS:1;orig_transcript_id=gnl|NXY56|LgyM4147_000005000.1;product=Protein of unknown function (DUF2946)%2C putative
# CP103914.1	Genbank	exon	3908	4906	.	-	.	ID=exon-LgyM4147_000005000.1-1;Parent=rna-LgyM4147_000005000.1;gbkey=mRNA;locus_tag=NXY56_000001;orig_protein_id=gnl|NXY56|LgyM4147_000005000.1:CDS:1;orig_transcript_id=gnl|NXY56|LgyM4147_000005000.1;product=Protein of unknown function (DUF2946)%2C putative
# CP103914.1	Genbank	CDS	3908	4906	.	-	0	ID=cds-XQJ24122.1;Parent=rna-LgyM4147_000005000.1;Dbxref=NCBI_GP:XQJ24122.1;Name=XQJ24122.1;gbkey=CDS;locus_tag=NXY56_000001;orig_transcript_id=gnl|NXY56|LgyM4147_000005000.1;product=Protein of unknown function (DUF2946)%2C putative;protein_id=XQJ24122.1

### tRNA
# CP103918.1	Genbank	gene	369209	369280	.	+	.	ID=gene-NXY56_000484;Name=NXY56_000484;gbkey=Gene;gene_biotype=tRNA;locus_tag=NXY56_000484
# CP103918.1	Genbank	tRNA	369209	369280	.	+	.	ID=rna-NXY56_000484;Parent=gene-NXY56_000484;gbkey=tRNA;locus_tag=NXY56_000484;product=tRNA-Arg
# CP103918.1	Genbank	exon	369209	369280	.	+	.	ID=exon-NXY56_000484-1;Parent=rna-NXY56_000484;gbkey=tRNA;locus_tag=NXY56_000484;product=tRNA-Arg

### rRNA
# CP103922.1	Genbank	gene	429912	430030	.	+	.	ID=gene-NXY56_001028;Name=NXY56_001028;gbkey=Gene;gene_biotype=rRNA;locus_tag=NXY56_001028
# CP103922.1	Genbank	rRNA	429912	430030	.	+	.	ID=rna-NXY56_001028;Parent=gene-NXY56_001028;gbkey=rRNA;locus_tag=NXY56_001028;product=5S ribosomal RNA
# CP103922.1	Genbank	exon	429912	430030	.	+	.	ID=exon-NXY56_001028-1;Parent=rna-NXY56_001028;gbkey=rRNA;locus_tag=NXY56_001028;product=5S ribosomal RNA

### ncRNA
# CP103915.1	Genbank	gene	267089	267116	.	+	.	ID=gene-NXY56_008321;Name=NXY56_008321;Note=slRNA;gbkey=Gene;gene_biotype=ncRNA;locus_tag=NXY56_008321
# CP103915.1	Genbank	ncRNA	267089	267116	.	+	.	ID=rna-NXY56_008321;Parent=gene-NXY56_008321;gbkey=ncRNA;locus_tag=NXY56_008321;product=spliced_leader_RNA
# CP103915.1	Genbank	exon	267089	267116	.	+	.	ID=exon-NXY56_008321-1;Parent=rna-NXY56_008321;gbkey=ncRNA;locus_tag=NXY56_008321;product=spliced_leader_RNA

## Example 2:
### CDS
# CP027540.1      Genbank gene    1420849 1421355 .       -       .       ID=gene-SPV_1401;Name=folA;Note=Corresponds to SPD_1401 found in INSD CP000410.1~Corresponds to spr1429 found in RefSeq NC_003098.1~Corresponds to SP_1571 found in NC_003028.3;gbkey=Gene;gene=folA;gene_biotype=protein_coding;gene_synonym=dfr,dhfR;locus_tag=SPV_1401
# CP027540.1      Genbank CDS     1420849 1421355 .       -       0       ID=cds-AVN86481.1;Parent=gene-SPV_1401;Dbxref=GO:0004146,GO:0006545,GO:0009165,GO:0050661,GO:0055114,InterPro:IPR001796,InterPro:IPR012259,InterPro:IPR017925,InterPro:IPR024072,PFAM:PF00186,NCBI_GP:AVN86481.1;Name=AVN86481.1;gbkey=CDS;gene=folA;locus_tag=SPV_1401;product=Dihydrofolate reductase;protein_id=AVN86481.1;transl_table=11

### tRNA
# CP027540.1      Genbank gene    14827   14898   .       +       .       ID=gene-SPV_0015;Name=tRNA-Glu-1;Note=Corresponds to SPD_0015 found in INSD CP000410.1~Corresponds to sprt01 found in RefSeq NC_003098.1~Corresponds to SP_2242 found in NC_003028.3;gbkey=Gene;gene=tRNA-Glu-1;gene_biotype=tRNA;locus_tag=SPV_0015
# CP027540.1      Genbank tRNA    14827   14898   .       +       .       ID=rna-SPV_0015;Parent=gene-SPV_0015;Dbxref=RFAM:RF00005;Note=tRNA-Glu-UUC~Supported by PubMed ID 27174935;gbkey=tRNA;gene=tRNA-Glu-1;locus_tag=SPV_0015;product=tRNA-Glu
# CP027540.1      Genbank exon    14827   14898   .       +       .       ID=exon-SPV_0015-1;Parent=rna-SPV_0015;Dbxref=RFAM:RF00005;Note=tRNA-Glu-UUC~Supported by PubMed ID 27174935;gbkey=tRNA;gene=tRNA-Glu-1;locus_tag=SPV_0015;product=tRNA-Glu

### rRNA
# CP027540.1      Genbank gene    15138   16699   .       +       .       ID=gene-SPV_0016;Name=rrsA;Note=Corresponds to SPD_0016 found in INSD CP000410.1~Corresponds to sprr01 found in RefSeq NC_003098.1~Corresponds to SP_rrnaA16S found in NC_003028.3;gbkey=Gene;gene=rrsA;gene_biotype=rRNA;locus_tag=SPV_0016
# CP027540.1      Genbank rRNA    15138   16699   .       +       .       ID=rna-SPV_0016;Parent=gene-SPV_0016;Dbxref=RFAM:RF00177;Note=Small subunit ribosomal RNA;gbkey=rRNA;gene=rrsA;locus_tag=SPV_0016;product=16S ribosomal RNA
# CP027540.1      Genbank exon    15138   16699   .       +       .       ID=exon-SPV_0016-1;Parent=rna-SPV_0016;Dbxref=RFAM:RF00177;Note=Small subunit ribosomal RNA;gbkey=rRNA;gene=rrsA;locus_tag=SPV_0016;product=16S ribosomal RNA

### ncRNA
# CP027540.1      Genbank gene    86084   86201   .       +       .       ID=gene-SPV_2577;Name=Spd_sr9;gbkey=Gene;gene=Spd_sr9;gene_biotype=ncRNA;locus_tag=SPV_2577
# CP027540.1      Genbank ncRNA   86084   86201   .       +       .       ID=rna-SPV_2577;Parent=gene-SPV_2577;gbkey=ncRNA;gene=Spd_sr9;locus_tag=SPV_2577;product=ncRNA of unknown function
# CP027540.1      Genbank exon    86084   86201   .       +       .       ID=exon-SPV_2577-1;Parent=rna-SPV_2577;gbkey=ncRNA;gene=Spd_sr9;locus_tag=SPV_2577;product=ncRNA of unknown function

def extract_feature_annotations(gff_file, feature_type, output_tsv):
    with open(output_tsv, "w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["gene_id", "chr", "start", "end", "strand", "size", "locus_tag", "gene_name", "gbkey", "description"])

        with open(gff_file) as gff:
            for line in gff:
                if line.startswith("#") or "\t" not in line:
                    continue

                fields = line.strip().split("\t")
                if len(fields) != 9:
                    continue
                chr, source, feature, start, end, score, strand, frame, attributes = fields

                start = int(fields[3])
                end = int(fields[4])
                size = end - start + 1

                attr_dict = parse_gff3_attributes(attributes)

                if feature in feature_type:
                    gene_id = attr_dict.get("ID")
                    locus_tag = attr_dict.get("locus_tag", "-")
                    gene_name = attr_dict.get("gene", "-")
                    gbkey = attr_dict.get("gbkey", "-")
                    description = attr_dict.get("product", "-")

                    writer.writerow([gene_id, chr, start, end, strand, size, locus_tag, gene_name, gbkey, description])

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="""
        Extract gene annotations from NCBI Genbank GFF3 file.
        Example usage: python gff_to_tsv_ncbi.py --gff_file input.gff3 --feature_types CDS tRNA rRNA ncRNA --output_tsv output.tsv
        """
        )
    
    parser.add_argument("--gff_file", required=True, help="NCBI Genbank GFF3 file")
    parser.add_argument("--feature_types", nargs="+", default=["gene"], help="List of feature type to extract (e.g., CDS tRNA rRNA ncRNA). Default is 'gene'")
    parser.add_argument("--output_tsv", required=True, help="Gene annotation output TSV file")

    args = parser.parse_args()

    extract_feature_annotations(args.gff_file, args.feature_types, args.output_tsv)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()


