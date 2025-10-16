#!/usr/bin/env python3
from Bio import SeqIO
import sys, os, glob

# --- Check arguments ---
if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_folder> <fasta_extension>")
    print(f"Example: {sys.argv[0]} fasta_files fasta")
    sys.exit(1)

folder = sys.argv[1]
extension = sys.argv[2].lstrip('.')  # remove leading dot if user adds .fna

# --- Header ---
print("Gene\tTotal_sequences\tUnique_nucleotide_alleles\tUnique_protein_alleles\tAA/NT_alleles_ratio")

# --- Loop through all FASTA files ---
for fasta_file in sorted(glob.glob(os.path.join(folder, f"*.{extension}"))):
    nt_seqs = set()
    aa_seqs = set()
    total = 0

    for record in SeqIO.parse(fasta_file, "fasta"):
        total += 1
        seq = str(record.seq).upper().replace("-", "").replace("N", "")
        if len(seq) % 3 != 0 or len(seq) < 3:
            # skip incomplete or frame-shifted sequences
            continue

        nt_seqs.add(seq)

        try:
            aa_seq = str(record.seq.translate(table=11, to_stop=True))
            aa_seqs.add(aa_seq)
        except Exception:
            continue

    gene = os.path.splitext(os.path.basename(fasta_file))[0]
    nt_count = len(nt_seqs)
    aa_count = len(aa_seqs)
    ratio = f"{aa_count/nt_count:.2f}" if nt_count > 0 else "NA"

    print(f"{gene}\t{total}\t{nt_count}\t{aa_count}\t{ratio}")
