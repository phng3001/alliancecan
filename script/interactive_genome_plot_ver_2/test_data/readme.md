# Merge tsv files
```python
python merge_tsv.py -i "Linfantum*" -c "#CHROM,START,END" -t "NORMALIZED_LOG2_READS_RATIO_MUTANT/WT" -o normalized_log2_read_ratio_mutant_vs_wt.tsv --sort
```

# Rename columns
```bash
sed -i 's/#CHROM/chr/' normalized_log2_read_ratio_mutant_vs_wt.tsv
sed -i 's/START/start/' normalized_log2_read_ratio_mutant_vs_wt.tsv
sed -i 's/END/end/' normalized_log2_read_ratio_mutant_vs_wt.tsv
```

# Make genome plot
```bash
apptainer run python_3.11.5.sif python genome_plot.py \
-v normalized_log2_read_ratio_mutant_vs_wt.tsv \
-f TriTrypDB-68_LinfantumJPCM5.gff \
-t protein_coding_gene \
--output genome_plot.html \
--feature-height 0.15 \
--value-ytitle "Normalized log2 read ratio mutant/WT" \
--feature-ytitle "CDS"
```
