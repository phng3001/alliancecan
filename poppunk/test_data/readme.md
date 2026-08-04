# Create `input_path_file`
```bash
for X in $(cat sample_list.txt); do echo -e "$X\tinput/${X}.fasta"; done > input_path.tab
```

# Launch script
```
sbatch run_apptainer_poppunk.sh poppunk_v2.7.8.sif input_path.tab output
```
