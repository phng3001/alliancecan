# File formats

## Sequence format

### FASTA
#### Purpose
Store DNA, RNA, or protein sequences

#### Common file extensions
* `.fasta`
* `.fa`
* `.fas`
* `.fna`: nucleotide sequences
* `.ffn`: nucleotide coding sequences (CDS)
* `.faa`: protein (amino acids) sequences


#### Description
Each record contains:
1. A header line beginning with `>` followed by the sequence identifier and optional description
2. One or more sequence lines, until the next `>` or end of file
> **Notes:**
* FASTA headers may contain structured fields, but there is no universal standard
* Sequences use IUPAC codes
* Sequence lines may be wrapped (e.g., 60/80 characters) or unwrapped (all on one line). Some bioinformatics tools require unwrapped sequences 

#### Example
```
>seq1 description
ATGCGTATAGCTAGCTAGTAA
>seq2 description
GTGCTANNNCTAGGCTAGTAG
```



### FASTQ
#### Purpose
Store sequencing reads with per-base quality scores, typically from sequencing instruments

#### Common file extensions
* `.fastq` or `.fastq.gz` (compressed)
* `.fq` or `.fq.gz` (compressed)

#### Description
Each record contains 4 lines:
1. Header starting with `@` (sequence identifier)
2. Sequence
3. Separator line starting with `+` (optional repeated identifier)
4. Quality string, ASCII-encoded
> **Notes:**
* Quality string must match sequence length (1 character per base)
* In paired-end sequencing FASTQ files are distributed in pairs: `R1` (forward) and `R2` (reverse)
* Quality scores represent the probability of sequencing error for each base:
```
Q = -10 × log10(P_error)
```
Quality score (Q) | Error probability (P_error) | Base call accuracy 
------------------|-----------------------------|-------------------
10 | 0.1 | 90%
20 | 0.01 | 99%
30 | 0.001 | 99.9%

* Nearly all modern FASTQ files (Illumina 1.8+) use Phred+33 (Sanger standard); Phred+64 is now obsolete
* Illumina CASAVA 1.8+ FASTQ header format:
```
@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<barcode(s)>
```

#### Example
```
@A01114:203:HGJKVDSXF:3:1101:1262:1016 1:N:0:GGCCAATATT+TATCGGACCG
GNGCATGACTAAACGTGTACTCATTACAGGAGTGAGTTCAGGGATCGGATTGGCTCAAGCTCGCCTCTTTTTAGAGAAGGGCTATCAAGTTTATGGAGTTGACCAAGGTGAAAAGCCACTCTTAGAGGGTGATTTTCGCTTTTTACAGAGA
+
F#FFFFFFFFFFFFFFFFFFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFF,:FFFFF:FF:,F:FFFF:FFFFFFFFFF
```
