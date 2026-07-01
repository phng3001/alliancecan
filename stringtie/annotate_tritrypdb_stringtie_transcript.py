#!/usr/bin/env python3

import argparse


# -----------------------------
# Attribute parsers
# -----------------------------

def parse_gtf_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(" ", 1)
        if len(parts) != 2:
            continue
        key, val = parts
        attrs[key] = val.replace('"', '').strip()
    return attrs


def parse_gff_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, val = item.split("=", 1)
            attrs[key] = val
    return attrs


# -----------------------------
# Load gene annotations
# -----------------------------

def load_gff_annotations(gff_file):

    annotations = {}
    genes = []

    with open(gff_file) as f:
        for line in f:

            if line.startswith("#") or "\t" not in line:
                continue

            cols = line.rstrip().split("\t")

            if cols[2] != "protein_coding_gene" and cols[2] != "ncRNA_gene":
                continue

            chrom = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            strand = cols[6]

            attrs = parse_gff_attributes(cols[8])

            locus_tag = attrs.get("ID", "")
            gene = attrs.get("Name", "-")

            key = (chrom, start, end, strand)

            annotations[key] = {
                "locus_tag": locus_tag,
                "gene": gene
            }

            genes.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "locus_tag": locus_tag,
                "gene": gene
            })

    return annotations, genes


# -----------------------------
# Build description dictionary
# -----------------------------

def build_description_dict(gff_file):

    desc_dict = {}

    with open(gff_file) as f:
        for line in f:

            if line.startswith("#") or "\t" not in line:
                continue

            if "ID=" not in line or "description=" not in line:
                continue

            cols = line.rstrip().split("\t")

            attrs = parse_gff_attributes(cols[8])

            locus_tag = attrs.get("ID")
            product = attrs.get("description")

            if locus_tag and product:
                desc_dict[locus_tag] = product

    return desc_dict


# -----------------------------
# Overlap detection
# -----------------------------

def find_overlapping_genes(chrom, start, end, strand, genes):

    overlaps = []

    for g in genes:

        if g["chr"] != chrom:
            continue

        if g["strand"] != strand:
            continue

        if not (end < g["start"] or start > g["end"]):
            overlaps.append(g)

    if not overlaps:
        return None

    locus_tags = " | ".join([g["locus_tag"] for g in overlaps])
    gene_names = " | ".join([g["gene"] if g["gene"] else "-" for g in overlaps])

    return locus_tags, gene_names


# -----------------------------
# Antisense detection
# -----------------------------

def find_antisense_genes(chrom, start, end, strand, genes):

    antisense_hits = []

    for g in genes:

        if g["chr"] != chrom:
            continue

        if g["strand"] == strand:
            continue

        if not (end < g["start"] or start > g["end"]):
            antisense_hits.append(g)

    if not antisense_hits:
        return None

    locus_tags = " | ".join([f"anti-{g['locus_tag']}" for g in antisense_hits])

    gene_names = " | ".join(
        [f"anti-{g['gene']}" if g["gene"] else "anti-" for g in antisense_hits]
    )

    return locus_tags, gene_names


# -----------------------------
# Flanking gene detection
# -----------------------------

def find_flanking_genes(chrom, start, end, genes):

    left_gene = None
    right_gene = None

    for g in genes:

        if g["chr"] != chrom:
            continue

        if g["end"] < start:
            if left_gene is None or g["end"] > left_gene["end"]:
                left_gene = g

        elif g["start"] > end:
            if right_gene is None or g["start"] < right_gene["start"]:
                right_gene = g

    if left_gene or right_gene:

        left_tag = left_gene["locus_tag"] if left_gene else "-"
        right_tag = right_gene["locus_tag"] if right_gene else "-"

        left_gene_name = left_gene["gene"] if left_gene and left_gene["gene"] else "-"
        right_gene_name = right_gene["gene"] if right_gene and right_gene["gene"] else "-"

        locus_tags = f"{left_tag} | {right_tag}"
        gene_names = f"{left_gene_name} | {right_gene_name}"

        return locus_tags, gene_names

    return None


# -----------------------------
# Description lookup
# -----------------------------

def get_description(locus_tag_string, desc_dict):

    if locus_tag_string in ["", "-"]:
        return "-"

    parts = [x.strip() for x in locus_tag_string.split("|")]

    descriptions = []

    for p in parts:

        is_antisense = False

        if p.startswith("anti-"):
            is_antisense = True
            p = p.replace("anti-", "")

        desc = desc_dict.get(p, "-")

        if is_antisense and desc != "-":
            desc = f"antisense {desc}"

        descriptions.append(desc)

    return " | ".join(descriptions)


# -----------------------------
# Transcript annotation
# -----------------------------

def annotate_transcripts(gtf_file, annotations, genes, desc_dict, output):

    with open(output, "w") as out:

        out.write("transcript_id\tchr\tstart\tend\tstrand\tlocus_tag\tgene\tdescription\n")

        with open(gtf_file) as f:

            for line in f:

                if line.startswith("#") or "\t" not in line:
                    continue

                cols = line.rstrip().split("\t")

                if cols[2] != "transcript":
                    continue

                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]

                attrs = parse_gtf_attributes(cols[8])
                transcript_id = attrs.get("transcript_id", "")

                key = (chrom, start, end, strand)

                # 1 exact coordinate match
                if key in annotations:

                    locus_tag = annotations[key]["locus_tag"]
                    gene = annotations[key]["gene"]

                else:

                    # 2 extended transcript
                    overlap = find_overlapping_genes(chrom, start, end, strand, genes)

                    if overlap:
                        locus_tag, gene = overlap

                    else:

                        # 3 antisense transcript
                        antisense = find_antisense_genes(chrom, start, end, strand, genes)

                        if antisense:
                            locus_tag, gene = antisense

                        else:

                            # 4 intergenic transcript
                            flank = find_flanking_genes(chrom, start, end, genes)

                            if flank:
                                locus_tag, gene = flank
                            else:
                                locus_tag = "-"
                                gene = "-"

                description = get_description(locus_tag, desc_dict)

                out.write(
                    f"{transcript_id}\t{chrom}\t{start}\t{end}\t{strand}\t{locus_tag}\t{gene}\t{description}\n"
                )


# -----------------------------
# Main
# -----------------------------

def main():

    parser = argparse.ArgumentParser(
        description="Annotate StringTie transcripts using TriTrypDB GFF"
    )

    parser.add_argument(
        "-g", "--gtf",
        required=True,
        help="StringTie merged GTF"
    )

    parser.add_argument(
        "-a", "--gff",
        required=True,
        help="TriTrypDB reference GFF"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV"
    )

    args = parser.parse_args()

    annotations, genes = load_gff_annotations(args.gff)

    desc_dict = build_description_dict(args.gff)

    annotate_transcripts(
        args.gtf,
        annotations,
        genes,
        desc_dict,
        args.output
    )


if __name__ == "__main__":
    main()

