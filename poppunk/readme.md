# Source
https://github.com/bacpop/PopPUNK

# Function
Alignment-free, variable-length k-mer comparisons to rapidly cluster genomes into distinct strains for epidemiological tracking

# Container
poppunk_v2.7.8.sif

# Usage
```bash
bash run_apptainer_poppunk.sh <container> <input_path_file> <output_dir>
```
- `input_path_file`: Tab-delimited file listing query input assembly path, one assembly per line. Cf. example in the `test_data` folder 
- `output_dir`: Name for the output directory. It will also be used as prefix for output files. For strain typing, refer to the output file with the suffix "_external_clusters.csv".
- This script was written for *S. pneumoniae* typing using the GPS database version 11. For other species or other database versions,  the `poppunk_data`, `distances` and `external_clustering` variables need to be edited.

# Notes
It seems that speed depends more on the database than the number of input strains: On Nibi compute nodes it tooks ~ 4 min for 3 *S. pneumoniae* strains and ~ 9 min for 1500 strains.

# More info
https://poppunk-docs.bacpop.org/
