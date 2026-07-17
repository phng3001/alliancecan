# Source
https://github.com/davidemms/OrthoFinder
https://github.com/OrthoFinder/OrthoFinder

# Container
orthofinder_v3.1.5.sif

# Usage
## 1. Prepare the input
- OrthoFinder requires **protein** sequences, not nucleotide sequences. Each genome should have one protein FASTA file.

E.g.
```
proteomes/
├── Lguyanensis.faa
├── Linfantum.faa
└── Lmajor.faa
```
> **Notes:** Each protein ID should be unique.

## 2. Run OrthoFinder 
```bash
bash run_apptainer_orthofinder.sh <container> <input_dir_path>
```
OrthoFinder output (`OrthoFinder`) will be written in the same directory as the input one containing the proteomes
