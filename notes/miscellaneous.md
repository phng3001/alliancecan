# Create a shared directory
```bash
mkdir shared_dir_name
chgrp def-mouellet shared_dir_name # def-professor
chmod -R 2775 shared_dir_name # owner = rwx, group = rwx, others = r-x
setfacl -R -m g::rwx shared_dir_name
setfacl -d -m g::rwx shared_dir_name
setfacl -d -m o::r-x shared_dir_name
```

# Bash script: disable block
```bash
: <<'END'
echo "This line will not be executed"
END
```

# Data stream
Operator | Function | Syntax
---------|----------|-------
`>` | Redirect stdout | `command > file.txt` or `command 1> file.txt`
`2>` | Redirect stderr | `command 2> error.txt`
`&>` | Redirect all | `command &> log.txt`
`>>` | Append stdout | `command >> files.txt`
`2>>` | Append stderr | `command 2>> errors.txt`
`&>>` | Append all | `command &>> log.txt`

> **Notes:** Redirect to `/dev/null` to hide the corresponding stream (E.g. `command &> /dev/null`)

# Task status marks
Symbol | Code | Meaning
-------|------|--------
 ✅ | "\u2705" | Success / Done
 ❌ | "\u274C" | Failure
 ⚠ | "\u26A0" | Warning
 ☐ | "\u2610" | Pending / Not done

Color | Code
------|-----
Green | "\e[32m"
Red | "\e[31m"
Yellow | "\e[33m"

E.g.
```bash
# ⚠ in yellow
echo -e "\e[33m\u26A0" 
```

# Find files not belonging to specific users
```bash
find /path/to/directory -type f ! -user user1 ! -user user2
```

# Find files belonging to specific group
```bash
find /path/to/directory -group somegroupname
```

# Change group ownership recursively
```bash
chown -h -R owner:groupname /path/to/directory
```

# Google Colab
## Upload
```python
from google.colab import files
upload = files.upload()
```

## Download
```python
from google.colab import files
download = files.download('output.csv')
```

# Fasta header formatting
```bash
# tritrypdb
# >LINF_010005000-T1-p1 | transcript=LINF_010005000-T1 | gene=LINF_010005000 | organism=Leishmania_infantum_JPCM5 | gene_product=Protein of unknown function (DUF2946) | transcript_product=Protein of unknown function (DUF2946) | location=LinJ.01:3710-4711(-) | protein_length=333 | sequence_SO=chromosome | SO=protein_coding_gene | is_pseudo=false
sed -E 's/^>.*gene=([^ ]+).*/>\1/' input.fasta > output.fasta
# >LINF_010005000
```

```bash
# genbank
# >lcl|CP103914.1_prot_XQJ24122.1_1 [locus_tag=NXY56_000001] [protein=Protein of unknown function (DUF2946), putative] [protein_id=XQJ24122.1] [location=complement(3908..4906)] [gbkey=CDS]
sed -E 's/^>.*\[locus_tag=([^]]+)\].*/>\1/' input.fasta > output.fasta
# >NXY56_000001
```

# Parsing GFF/GTF
## GFF
```python
def parse_gff3_attributes(attr_string):
    """Parse GFF attribute string (9th column of a GFF file) into a dictionary."""
    attr_dict = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attr_dict[key] = value
    return attr_dict
```

## GTF
```python
def parse_gtf_attributes(attr_string):
    """Parse GTF attribute string (9th column of a GTF file) into a dictionary."""
    attr_dict = {}
    for key, value in re.findall(r'\s*([^ ;]+)\s+"([^"]+)"', attr_string):
        attrs[key] = value
    return attr_dict
```

# Replace string in each filename
```bash
# Replace 'genes' with 'features' in filenames
for f in *genes*; do mv -- "$f" "${f//genes/features}"; done
```
> Preview before renaming:
```bash
for f in *genes*; do echo mv -- "$f" "${f//genes/features}"; done
```

# Remove the first N lines
## tail
```bash
tail -n +$((N+1)) file.txt
```
## sed
```bash
sed '1,Nd' file.txt
```
## awk
```bash
awk 'NR>N' file.txt
```

# Merge files with the same header
```bash
awk 'FNR==1{if(!hdr){hdr=$0; print; next} if($0!=hdr){print "Headers differ!" > "/dev/stderr"; exit}} FNR>1' *.tsv > merged.tsv
```

# seqkit
```bash
# Filter out contigs < 500bp
seqkit seq -m 500 input.fasta > output.fasta
```
