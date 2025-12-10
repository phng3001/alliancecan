# File formats

## Sequence format

### FASTA
#### Purpose
Store DNA, RNA, or protein sequences

#### Description
##### Each record contains:
* A header line beginning with `>` followed by the sequence identifier and optional description
* One or more sequence lines, until the next `>` or end of file
##### Notes:
* FASTA headers may contain structured fields, but there is no universal standard
* Sequences use IUPAC codes
* Sequence lines may be wrapped (e.g., 60/80 characters) or unwrapped (all on one line). Some bioinformatics tools require unwrapped sequences 

#### Example
```bash
>seq1 description
ATGCGTATAGCTAGCTAGTAA
>seq2 description
GTGCTANNNCTAGGCTAGTAG
```