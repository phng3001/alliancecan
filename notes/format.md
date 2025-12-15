---
title:  "File formats"
excerpt: "Common file formats used in bioinformatics"
collection: cheatsheets
---

- [Sequence format](#sequence-format)
  - [FASTA](#fasta)
    - [Purpose](#purpose)
    - [Common file extensions](#common-file-extensions)
    - [Description](#description)
    - [Example](#example)
  - [FASTQ](#fastq)
    - [Purpose](#purpose-1)
    - [Common file extensions](#common-file-extensions-1)
    - [Description](#description-1)
    - [Example](#example-1)
- [Annotation format](#annotation-format)
  - [BED (Browser Extensible Data)](#bed-browser-extensible-data)
    - [Purpose](#purpose-2)
    - [Common file extensions](#common-file-extensions-2)
    - [Description](#description-2)
    - [Example](#example-2)
      - [BED3](#bed3)
      - [BED6](#bed6)
      - [BED12](#bed12)
  - [GFF3 (General Feature Format)](#gff3-general-feature-format)
    - [Purpose](#purpose-3)
    - [Common file extensions](#common-file-extensions-3)
    - [Description](#description-3)
    - [Example](#example-3)
  - [GTF (Gene Transfer Format)](#gtf-gene-transfer-format)
    - [Purpose](#purpose-4)
    - [Common file extensions](#common-file-extensions-4)
    - [Description](#description-4)
    - [Example](#example-4)
- [Alignment format](#alignment-format)
  - [SAM / BAM / CRAM](#sam--bam--cram)
    - [Purpose](#purpose-5)
    - [Common file extensions](#common-file-extensions-5)
    - [Description](#description-5)
      - [The header section](#the-header-section)
      - [The alignment section](#the-alignment-section)
    - [Example](#example-5)
  - [VCF / BCF](#vcf--bcf)
    - [Purpose](#purpose-6)
    - [Common file extensions](#common-file-extensions-6)
    - [Description](#description-6)
      - [Meta-information lines](#meta-information-lines)
      - [Header line](#header-line)
      - [Data lines](#data-lines)
        - [`INFO` (Column 8)](#info-column-8)
        - [`FORMAT` (Column 9)](#format-column-9)
    - [Example](#example-6)
- [More info](#more-info)


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
* In paired-end sequencing, FASTQ files are distributed in pairs: `R1` (forward) and `R2` (reverse)
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

**BED columns**

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
* Feature lines: exactly 9 **tab-delimited** columns
* Optional FASTA section starts with `##FASTA`
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
3 | type | Feature type (typically Sequence Ontology terms, e.g. `gene`, `exon`, `CDS`, `tRNA`, `ncRNA`)
4 | start | Start position (1-based, inclusive)
5 | end | End position (1-based, inclusive)
6 | score | Numeric score or `.` (no score)
7 | strand | `+`, `-` or `.` (no strand)
8 | phase | For CDS: `0`/`1`/`2` (reading frame). For non-CDS: `.`
9 | attributes | Semicolon-separated `key=value` pairs (e.g. `ID=gene1;Parent=transcript1`)

> **Notes:**
* Coordinates are **1-based, fully closed**: both `start` and `end` are included
* Must be **tab-delimited**, all 9 columns required
* Supports **multi-level hierarchies** using `ID=` and `Parent=` attributes. E.g. gene → mRNA → exon/CDS
* Attribute values must be URL-encoded if they contain special characters (e.g. spaces → `%20`, `;` → `%3B`, `=` → `%3D`)
* Different databases and annotation pipelines (e.g. NCBI, Ensembl, Prokka) also have their **own conventions**, e.g. in feature types, attribute keys

#### Example
```
##gff-version 3
chr1    Ensembl    gene        1000    5000    .    +    .    ID=gene1;Name=GeneA
chr1    Ensembl    mRNA        1000    5000    .    +    .    ID=transcript1;Parent=gene1
chr1    Ensembl    exon        1000    1200    .    +    .    ID=exon1;Parent=transcript1
chr1    Ensembl    exon        2000    2200    .    +    .    ID=exon2;Parent=transcript1
chr1    Ensembl    CDS         1050    1200    .    +    0    ID=cds1;Parent=transcript1
chr1    Ensembl    CDS         2000    2150    .    +    2    ID=cds2;Parent=transcript1
```



### GTF (Gene Transfer Format)
#### Purpose
Describe gene structure and annotations (e.g. genes, transcripts, exons, UTRs, CDS), primarily for gene-centric analysis like RNA-seq and gene quantification

#### Common file extensions
* `.gtf`: recommended standard
* `.gff`: sometimes used, but ambiguous (may also be GFF3)

#### Description
A GTF file is very similar to GFF3 but uses a **different attribute syntax** (column 9) and follows a **GFF2-derived specification**

A GTF file includes:
* Optional header lines (start with `#`)
* Feature lines: exactly 9 **tab-delimited** columns

**GTF columns**

Column | Name | Description
-------|------|------------
1 | seqname | Sequence name (chromosome/contig ID)
2 | source | Annotation source (e.g. `RefSeq`, `Ensembl`, `StringTie`)
3 | feature | Feature type (e.g. `gene`, `transcript`, `exon`, `CDS`, `UTR`)
4 | start | Start position (1-based, inclusive)
5 | end | End position (1-based, inclusive)
6 | score | Numeric score or `.` (no score)
7 | strand | `+`, `-` or `.` (no strand)
8 | frame | For CDS: `0`/`1`/`2` (reading frame). For non-CDS: `.`
9 | attributes | Semicolon-separated `key "value"` pairs (e.g. `gene_id "G1"; transcript_id "T1"`)

> **Notes:**
* Coordinates are **1-based, fully closed**: both `start` and `end` are included
* Must be **tab-delimited**, all 9 columns required
* GTF does **not support embedded FASTA** section
* Hierarchies are **implicit** using `gene_id` and `transcript_id` attributes
* Attribute values are **quoted** strings
* Common attributes: 
    - `gene_id` (required)
    - `transcript_id` (required for transcript-level features)
    - `gene_name`, `gene_type`, `gene_biotype`, `transcript_name`, `exon_number` (optional)
* Most RNA-seq tools expect GTF, not GFF3

#### Example
```
chr1    Ensembl    gene        1000    5000    .    +    .    gene_id "gene1"; gene_name "GeneA";
chr1    Ensembl    transcript  1000    5000    .    +    .    gene_id "gene1"; transcript_id "transcript1";
chr1    Ensembl    exon        1000    1200    .    +    .    gene_id "gene1"; transcript_id "transcript1"; exon_number "1";
chr1    Ensembl    exon        2000    2200    .    +    .    gene_id "gene1"; transcript_id "transcript1"; exon_number "2";
chr1    Ensembl    CDS         1050    1200    .    +    0    gene_id "gene1"; transcript_id "transcript1";
chr1    Ensembl    CDS         2000    2150    .    +    2    gene_id "gene1"; transcript_id "transcript1";
```



## Alignment format 
### SAM / BAM / CRAM
#### Purpose
Store read-to-reference alignments produced by sequence aligners, including alignment position, mapping quality, CIGAR operations, and optional tags

#### Common file extensions
* `.sam`: plain text, human-readable
* `.bam`: binary, compressed version of SAM
* `.cram`: reference-based, highly compressed format (smaller than BAM)

#### Description
SAM (Sequence Alignment Map), BAM (Binary Alignment Map), and CRAM (Compressed Reference-oriented Alignment Map) represent the same alignment information in different encodings: SAM is text-based, while BAM and CRAM are binary formats optimized for storage and performance

A SAM file includes:
* Optional header lines (start with `@`)
* Alignment records: one per read, with 11 mandatory columns plus optional fields

SAM is **tab-delimited** and designed to represent how sequencing reads align to **a reference genome**

##### The header section
* Each line begins with `@` followed by a two-letter header record type
* Each line is tab-delimited
* Except `@CO` lines, each data field follows a `TAG:VALUE` format

**Common record types and tags**

Record type | Description | Tags
------------|-------------|-----
`@HD` | File-level metadata (e.g. format version, sorting order) | `VN` (format), `SO` (sorting order)
`@SQ` | Reference sequence dictionary | `SN` (reference sequence name), `LN` (reference sequence length)
`@RG` | Read group | `ID` (read group identifier), `PL` (platform/technology), `SM` (sample)
`@PG` | Program information | `ID` (program record identifier), `PN` (program name), `VN` (program version)
`@CO` | Free-text comment

##### The alignment section
* Each line typically represents the linear alignment of a read to the reference genome
* Each line consists of 11 or more tab-delimited fields
    - The first 11 fields are mandatory
    - Optional fields start at column 12 and follow the format `TAG:TYPE:VALUE`

**SAM columns**

Column | Name | Description
-------|------|------------
1 | QNAME | Query/read name
2 | FLAG | Bitwise flag describing read properties
3 | RNAME | Reference sequence name
4 | POS | 1-based leftmost alignment position
5 | MAPQ | Mapping quality (Phred-scaled)
6 | CIGAR | Alignment operations (e.g. matches, insertions, deletions)
7 | RNEXT | Mate reference name (`=` if same as RNAME)
8 | PNEXT | Mate alignment position
9 | TLEN | Template length
10 | SEQ | Read sequence (`*` if not stored)
11 | QUAL | Base quality string (same as the quality string in FASTQ format)

> **Notes:**
* All lines in SAM are **tab-delimited**, including header lines
* Coordinates are **1-based**
* A read may occupy multiple alignment lines (secondary or supplementary alignments)
* Unmapped reads have `RNAME` set as `*` and `POS` as `0`
* BAM and CRAM require an **index** (`.bai`, `.csi`, `.crai`) for random access
* CRAM files require access to the reference genome used for compression
* Some tools require **coordinate-sorted** BAM/CRAM

#### Example
```
@HD	VN:1.6	SO:coordinate
@SQ	SN:ref	LN:45
r001	99	ref	7	30	8M2I4M1D3M	=	37	39	TTAGATAAAGGATACTG	*
r002	0	ref	9	30	3S6M1P1I4M	*	0	0	AAAAGATAAGGATA	*
r001	147	ref	37	30	9M	=	7	-39	CAGCGGCAT	*	NM:i:1
```


### VCF / BCF
#### Purpose
Store genomic variants (SNPs, indels, and structural variants) relative to a reference genome, together with quality metrics, filters, annotations, and per-sample genotypes

#### Common file extensions
* `.vcf`: plain text, human-readable
* `.vcf.gz`: bgzip-compressed VCF
* `.bcf`: binary version of VCF

#### Description
VCF (Variant Call Format) represents genomic variation at specific genomic positions relative to **a reference genome**. BCF (Binary Variant Call Format) stores the same information in a compact binary form

A VCF file includes:
* Meta-information lines (start with `##`)
* A header line: starts with `#CHROM`, **tab-delimited**
* Variant records: one variant per line, **tab-delimited**

##### Meta-information lines
* Start with `##`
* May be unstructured `##key=value` or structured `##key=<key1=value1,key2=value2,key3=value3>`

##### Header line
* Is tab-delimited
* Contains 8 mandatory columns:
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
```
* If genotype data is present, the 9th column is `FORMAT`, followed by one or more sample IDs (duplicate sample IDs are not allowed)

##### Data lines
* Are tab-delimited
* Missing values are specified with a `.`

**VCF columns**

Column | Name | Description
-------|------|------------
1 | CHROM | Chromosome/contig reference sequence name
2 | POS | 1-based variant position
3 | ID | Variant identifier
4 | REF | Reference allele
5 | ALT | Alternate allele(s), comma-separated
6 | QUAL | Variant quality score (Phred-scaled)
7 | FILTER | Filter status (`PASS`, or a semicolon-separated list of filters that fail)
8 | INFO | Semicolon-separated `key=value` annotations
9 | FORMAT | Colon-separated genotype field keys
10+ | Samples | Per-sample genotype data, following the `FORMAT` layout

###### `INFO` (Column 8)
* Describes **variant-level (site-level)** information shared by all samples
* Common keys:
    - `DP`: Total read depth across samples
    - `AF`: Allele frequency
    - `AC`: Allele count
    - `AN`: Total number of alleles
    - `MQ`: Mapping quality
    - `ANN`: Functional annotation (gene, effect)

###### `FORMAT` (Column 9)
* Describes the structure of per-sample fields
* Combines with sample column(s) (Column 10+) to store **sample-level** information
* Common keys:
    - `GT`: Genotype
    - `DP`: Read depth
    - `AD`: Allelic depths
    - `GQ`: Genotype quality
    - `PL`: Phred-scaled genotype likelihoods
* `GT`:
    - Encodes the alleles carried by a sample at a variant site
    - Interpreted together with the `REF` and `ALT` columns
    - Alleles are encoded as integers: `0` → reference allele (`REF`), `1` → first alternate allele, `2` → second alternate allele, etc
    - Alleles are separated by `/` (unphased) or `|` (phased)
    - The number of alleles reflects sample ploidy (e.g. `0` → haploid, `0/1` → diploid, `0/0/1` → triploid)

**Common diploid genotypes**

GT | Interpretation
---|---------------
`0/0` | Homozygous reference
`0/1` | Heterozygous (REF / ALT1)
`1/1` | Homozygous alternate (ALT1 / ALT1)
`1/2` | Heterozygous with two different alternate alleles (ALT1 / ALT2)
`./.` | Missing genotype: both alleles unknown
`0/.` | Partially missing genotype: One REF allele, one unknown

> **Notes:**
* Coordinates are **1-based**
* File must be **tab-delimited**
* **Multiallelic variants** are represented in a **single record** (e.g. `ALT=A,C`)
* `INFO` and `FORMAT` fields must be declared in meta-information lines
* `INFO` fields are variant-centric, while `FORMAT` fields are sample-centric
* **Structural variants** may use **symbolic alleles** (e.g. `<DEL>`, `<INS>`, `<DUP>`, `<INV>`, `<CNV>`)
* Indexed `.vcf.gz` files (`.tbi` or `.csi`) enable random access; only **bgzip** compression (not gzip) is compatible with indexing
* BCF files typically use `.csi` indexes for random access

#### Example
```
##fileformat=VCFv4.3
##contig=<ID=chr1,length=248956422>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample1
chr1	1050	rs123	A	G	60	PASS	DP=35	GT	0/1
chr1	2000	.	CT	C	50	PASS	DP=20	GT	1/1
chr1	5000	.	N	<DEL>	40	PASS	SVTYPE=DEL;END=6200;SVLEN=-1200	GT	0/1
```



## More info
https://genome.ucsc.edu/FAQ/FAQformat.html

https://samtools.github.io/hts-specs/SAMv1.pdf

https://samtools.github.io/hts-specs/VCFv4.3.pdf