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

# CDS lines in TriTrypDB GFF3 files
# LinJ.01	VEuPathDB	protein_coding_gene	9038	11059	.	-	.	ID=LINF_010005200;Name=KIN13-1;description=Kinesin-13;ebi_biotype=protein_coding
# LinJ.01	VEuPathDB	mRNA	9038	11059	.	-	.	ID=LINF_010005200-T1;Parent=LINF_010005200;description=Kinesin-13;gene_ebi_biotype=protein_coding
# LinJ.01	VEuPathDB	exon	9038	11059	.	-	.	ID=exon_LINF_010005200-T1-E1;Parent=LINF_010005200-T1;gene_id=LINF_010005200
# LinJ.01	VEuPathDB	CDS	9038	11059	.	-	0	ID=LINF_010005200-T1-p1-CDS1;Parent=LINF_010005200-T1;gene_id=LINF_010005200;protein_source_id=LINF_010005200-T1-p1

# Non coding RNA lines in TriTrypDB GFF3 files
# tRNA
# LinJ.02	VEuPathDB	ncRNA_gene	235340	235431	.	+	.	ID=LINF_020010500;description=anticodon:gacgac|aaVal;ebi_biotype=tRNA
# LinJ.02	VEuPathDB	tRNA	235340	235431	.	+	.	ID=LINF_020010500-T1;Parent=LINF_020010500;description=anticodon:gacgac|aaVal;gene_ebi_biotype=tRNA
# LinJ.02	VEuPathDB	exon	235340	235431	.	+	.	ID=exon_LINF_020010500-T1-E1;Parent=LINF_020010500-T1;gene_id=LINF_020010500

# rRNA
# LinJ.05	VEuPathDB	ncRNA_gene	41335	43371	.	+	.	ID=LINF_050006300;description=nucleolar RNA helicase II - putative;ebi_biotype=rRNA
# LinJ.05	VEuPathDB	rRNA	41335	43371	.	+	.	ID=LINF_050006300-T1;Parent=LINF_050006300;description=nucleolar RNA helicase II - putative;gene_ebi_biotype=rRNA
# LinJ.05	VEuPathDB	exon	41335	43371	.	+	.	ID=exon_LINF_050006300-T1-E1;Parent=LINF_050006300-T1;gene_id=LINF_050006300

# ncRNA
# LinJ.05	VEuPathDB	ncRNA_gene	452835	452907	.	+	.	ID=LINF_050017600;description=snoRNA;ebi_biotype=misc_RNA
# LinJ.05	VEuPathDB	ncRNA	452835	452907	.	+	.	ID=LINF_050017600-T1;Parent=LINF_050017600;description=snoRNA;gene_ebi_biotype=misc_RNA
# LinJ.05	VEuPathDB	exon	452835	452907	.	+	.	ID=exon_LINF_050017600-T1-E1;Parent=LINF_050017600-T1;gene_id=LINF_050017600

def extract_gene_annotations(gff_file, output_tsv):
    with open(output_tsv, "w", newline="") as output:
        writer = csv.writer(output, delimiter="\t")
        writer.writerow(["gene_id", "chr", "start", "end", "strand", "size", "locus_tag", "gene_name", "biotype", "description"])

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

                if feature in ["protein_coding_gene", "ncRNA_gene"]:
                    gene_id = attr_dict.get("ID")
                    locus_tag = gene_id
                    gene_name = attr_dict.get("Name", "-")
                    biotype = attr_dict.get("ebi_biotype", "-")
                    description = attr_dict.get("description", "-")

                    writer.writerow([gene_id, chr, start, end, strand, size, locus_tag, gene_name, biotype, description])

def main():
    parser = argparse.ArgumentParser(description="Extract gene annotations from TriTrypDB GFF3 file.")
    parser.add_argument("--gff_file", required=True, help="TriTrypDB GFF3 file")
    parser.add_argument("--output_tsv", required=True, help="Gene annotation output TSV file")
    args = parser.parse_args()

    extract_gene_annotations(args.gff_file, args.output_tsv)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()


