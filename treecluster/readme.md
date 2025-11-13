# Source
https://github.com/niemasd/TreeCluster

# Function
TreeCluster infers clusters from phylogenetic tree

# Container
treecluster_v1.0.4.sif

# Usage
```bash
apptainer run treecluster_v1.0.4.sif TreeCluster.py -h
```

# Example
```bash
apptainer run treecluster_v1.0.4.sif \
TreeCluster.py \
-i tree.nwk \
-m max \
-t 0.01 \
-o clusters.tsv
```
