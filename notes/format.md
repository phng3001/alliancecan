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
3. Separator line starting with `+` (may optionally repeat the identifier)
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



## Annotation format 

### BED (Browser Extensible Data)
#### Purpose
Define genomic regions of interest as coordinates (e.g. features, coverage windows, capture regions)

#### Common file extensions
* `.bed`

#### Description
A BED file is a **tab-separated** text file describing genomic intervals
* Column 1-3 are mandatory (BED3)
* Up to 12 columns are defined (BED12, commonly used for exon/intron structure)

Column | Name | Description
-------|------|------------
1 | chrom | Chromosome name
2 | chromStart | Start position (0-based, inclusive)
3 | chromEnd | End position (0-based, exclusive)
4 | name | Feature name
5 | score | Integer 0–1000
6 | strand | `+`, `-` or `.` (no strand)
7 | thickStart | Start position of thick drawing (browser display)
8 | thickEnd | End position of thick drawing (browser display)
9 | itemRgb | RGB color (e.g. `255,0,0`) (browser display)
10 | blockCount | Number of blocks (exons)
11 | blockSizes | Comma-separated list of block sizes. The number of items in this list should correspond to `blockCount`
12 | blockStarts | Comma-separated list of block starts relative to `chromStart`. The number of items in this list should correspond to `blockCount`

> **Notes:**
* Coordinates are **0-based, half-open**: `chromStart` is included, `chromEnd` is excluded
* Must be **tab-delimited**; columns cannot be skipped (i.e. BED6 requires columns 1–6, not 1–5 + 7)
* BED has **no required header**. Optional UCSC `track`/`browser` lines or comment lines `#` may appear
* Many tools expect **sorted** BED files (by `chrom`, then by `chromStart`)

#### Example
##### BED3
```
chr1    1000    2000
```
##### BED6
```
chr1    1000    2000    transcript1    900    +
```
##### BED12
```
chr1    1000    2000    transcript1    900    +    1000    2000    0,128,255    2    200,300    0,700
```



### GFF3 (General Feature Format)
#### Purpose
Describe genomic features and annotations (genes, transcripts, exons, CDS, ncRNAs, repeats, regulatory regions, etc)

#### Common file extensions
* `.gff3`: recommended standard
* `.gff`: sometimes GFF3, sometimes older GFF2

#### Description
A GFF3 file includes:
* Optional directives (header lines start with `##`, e.g. `##gff-version 3`, `##sequence-region`)
* Features lines: exactly 9 tab-delimited columns
* Optional FASTA section start with `##FASTA`
```
##FASTA
>scaffold1.1
ACAGCGCTGTCTTTGTGGCCTGCAGGG
```

**GFF3 columns**
Column | Name | Description
-------|------|------------
1 | seqid | Sequence name (chromosome/contig ID)
2 | source | Annotation source (e.g. `RefSeq`, `Ensembl`, `Prokka`)
3 | type | Feature type (e.g. `gene`, `exon`, `CDS`, `tRNA`, `ncRNA`)
4 | start | Start position (1-based, inclusive)
5 | end | End position (1-based, inclusive)
6 | score | Numeric score or `.` (no score)
7 | strand | `+`, `-` or `.` (no strand)
8 | phase | For CDS: `0`/`1`/`2` (reading frame). For non-CDS: `.`
9 | attributes | Semicolon-separated `key=value` pairs (e.g. `ID=gene1;Parent=transcript1`)

> **Notes:**
* Coordinates are `1-based, fully closed`: both `start` and `end` are included
* Must be `tab-delimited`, all 9 columns required
* Supports multi-level hierarchies using `ID=` and `Parent=` attributes. E.g. gene → mRNA → exon/CDS

#### Example
```
##gff-version 3
chr1    Ensembl    gene        1000    5000    .    +    .    ID=gene1;Name=GeneA
chr1    Ensembl    mRNA        1000    5000    .    +    .    ID=transcript1;Parent=gene1
chr1    Ensembl    exon        1000    1200    .    +    .    Parent=transcript1
chr1    Ensembl    exon        2000    2200    .    +    .    Parent=transcript1
chr1    Ensembl    CDS         1050    1200    .    +    0    Parent=transcript1
chr1    Ensembl    CDS         2000    2150    .    +    2    Parent=transcript1
```



## More info
https://genome.ucsc.edu/FAQ/FAQformat.html