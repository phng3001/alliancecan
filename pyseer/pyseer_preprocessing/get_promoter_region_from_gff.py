import sys
import argparse
import re

def parse_gff3_attributes(attr_string):
    """Parse GFF attribute string (9th column of a GFF file) into a dictionary."""
    attr_dict = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attr_dict[key] = value
    return attr_dict

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Create a BED file containing the intergenic region "
            "upstream of each CDS, extending from the CDS to the "
            "nearest upstream CDS."
        )
    )

    parser.add_argument("--gff", required=True, help="Input GFF3 file")
    parser.add_argument("--output", required=True, help="Output BED file")
    parser.add_argument("--max-upstream", type=int, default=None,
        help=(
            "Optional maximum promoter length. "
            "If omitted, promoter extends all the way to the "
            "nearest upstream CDS."
        )
    )

    args = parser.parse_args()

    # ============================================================
    # 1. Read all CDS features
    # ============================================================

    cds_features = []

    with open(args.gff) as gff:
        for line in gff:
            if line.startswith("#") or "\t" not in line:
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            chr, source, feature, start, end, score, strand, frame, attributes = fields

            # Get CDS features
            if feature != "CDS":
                continue

            start = int(start)
            end = int(end)

            attr_dict = parse_gff3_attributes(attributes)

            # gene_id = attr_dict.get("ID")
            locus_tag = attr_dict.get("locus_tag", "-")
            gene_name = attr_dict.get("gene", "-")

            cds_features.append({
                "chrom": chr,
                "start": start,
                "end": end,
                "strand": strand,
                "locus_tag": locus_tag,
                "gene_name": gene_name,
                "attributes": attributes
            })

    print(f"CDS features found: {len(cds_features)}")

    # ============================================================
    # 2. Generate upstream intergenic regions
    # ============================================================

    n_written = 0
    n_zero = 0
    n_clipped = 0

    with open(args.output, "w") as outfile:
        for cds in cds_features:
            chrom = cds["chrom"]
            start = cds["start"]
            end = cds["end"]
            strand = cds["strand"]
            locus_tag = cds["locus_tag"]

            # Find CDSs that are upstream or overlap the upstream side of the current CDS.

            # ----------------------------------------------------
            # PLUS STRAND
            # ----------------------------------------------------

            if strand == "+":

                # Upstream is toward smaller coordinates.

                upstream_cds = [
                    other
                    for other in cds_features
                    if (
                        other["chrom"] == chrom
                        and other is not cds
                        and other["start"] < start
                    )
                ]

                if upstream_cds:

                    # Nearest upstream CDS = the one with the largest end coordinate among upstream CDSs.
                    nearest = max(
                        upstream_cds,
                        key=lambda x: x["end"]
                    )

                    # Intergenic region starts immediately after the upstream CDS.
                    promoter_start = nearest["end"] + 1

                    # Intergenic region ends immediately before current CDS.
                    promoter_end = start - 1
                    clipped_by = nearest["locus_tag"]

                else:

                    # No upstream CDS on this chromosome.
                    # Extend to the beginning of the chromosome.
                    promoter_start = 1
                    promoter_end = start - 1
                    clipped_by = "chromosome_start"

            # ----------------------------------------------------
            # MINUS STRAND
            # ----------------------------------------------------

            elif strand == "-":

                # Upstream is toward larger coordinates.

                upstream_cds = [
                    other
                    for other in cds_features
                    if (
                        other["chrom"] == chrom
                        and other is not cds
                        and other["end"] > end
                    )
                ]

                if upstream_cds:

                    # Nearest upstream CDS = the one with the smallest end coordinate among upstream CDSs.

                    nearest = min(
                        upstream_cds,
                        key=lambda x: x["end"]
                    )

                    # Intergenic region starts immediately after current CDS.
                    promoter_start = end + 1

                    # Intergenic region ends immediately before upstream CDS.
                    promoter_end = nearest["start"] - 1
                    clipped_by = nearest["locus_tag"]

                else:

                    # No upstream CDS on this chromosome.
                    # Extend to the end of the chromosome.
                    promoter_start = end + 1
                    promoter_end = None
                    clipped_by = "chromosome_end"

            else:

                print(
                    f"WARNING: skipping CDS {locus_tag} "
                    f"with unknown strand '{strand}' "
                    f"at {chrom}:{start}-{end}"
                )

                continue

            # ====================================================
            # 3. Apply optional maximum upstream distance
            # ====================================================

            if args.max_upstream is not None:

                if strand == "+":

                    # Keep at most N bp upstream of start.
                    maximum_start = max(
                        1,
                        start - args.max_upstream
                    )

                    promoter_start = max(
                        promoter_start,
                        maximum_start
                    )

                elif strand == "-":

                    maximum_end = end + args.max_upstream

                    if promoter_end is None:
                        promoter_end = maximum_end
                    else:
                        promoter_end = min(
                            promoter_end,
                            maximum_end
                        )

            # ====================================================
            # 4. Handle zero-length promoters
            # ====================================================

            if promoter_end is not None:

                if promoter_start > promoter_end:

                    promoter_length = 0
                    n_zero += 1

                    continue

                promoter_length = (
                    promoter_end - promoter_start + 1
                )

            else:

                # No upstream CDS and no maximum distance.
                # We cannot create a valid BED interval without knowing chromosome length.
                # Therefore skip this promoter and report it.

                print(
                    f"WARNING: no upstream CDS found for "
                    f"{locus_tag} on {chrom}; "
                    f"chromosome length is required to extend "
                    f"to chromosome end."
                )

                continue

            # ====================================================
            # 5. Track whether region was limited
            # ====================================================

            if args.max_upstream is not None:

                if promoter_length < args.max_upstream:
                    n_clipped += 1

            # ====================================================
            # 6. Convert GFF coordinates to BED
            # ============================================================

            # GFF:
            #   1-based inclusive
            #
            # BED:
            #   0-based half-open

            bed_start = promoter_start - 1
            bed_end = promoter_end

            name = f"{locus_tag}_promoter"

            # ====================================================
            # 7. Write BED
            # ====================================================

            outfile.write(
                f"{chrom}\t"
                f"{bed_start}\t"
                f"{bed_end}\t"
                f"{name}\t"
                f"0\t"
                f"{strand}\t"
                f"{locus_tag}\t"
                f"{cds['gene_name']}\t"
                f"{promoter_length}\t"
                f"{clipped_by}\n"
            )

            n_written += 1

    # ============================================================
    # 8. Summary
    # ============================================================

    print()
    print("Summary")
    print("-------")
    print(f"CDS features found:         {len(cds_features)}")
    print(f"Promoters written:          {n_written}")
    print(f"Zero-length promoters:      {n_zero}")
    print(f"Promoters < maximum length: {n_clipped}")
    print(f"Output:                     {args.output}")



if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append('-h')
    main()
