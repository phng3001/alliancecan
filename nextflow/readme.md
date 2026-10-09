# Nextflow on Allicancecan's cluster
## Load Nextflow
```bash
module load StdEnv/2023 nextflow/26.04.4
```

## Setup Nextflow environment variables
|Variable|Description|Recommendation|
|--------|-----------|--------------|  
`NXF_HOME`|Nextflow's home directory, including downloaded framework files, plugins, and pipeline assets|A persistent location in your home or project space
`NXF_WORK`|Nextflow's work directory|Your scratch space if the number of files in project (shared among your group) is limiting. Otherwise use -work-dir or workDir in configuration instead
`NXF_SINGULARITY_CACHEDIR`|Cache directory for Singularity/Apptainer container images when using the Singularity integration|A shared, persistent storage of your group if you share pipelines. Otherwise a persistent storage in your project space, as the home storage quota in Alliancecan's cluster is limiting (only 50 GB)
`NXF_APPTAINER_CACHEDIR`|Cache directory for Apptainer images when using the Apptainer integration|Same as `NXF_SINGULARITY_CACHEDIR`

## Nextflow configuration
* Nextflow configuration can be managed using a Nextflow configuration file
* Nextflow configuration files are plain text files containing a set of properties
* The priority of a setting (highest to lowest):
    * Parameters specified on the command line (--param_name value)
    * Parameters provided using the -params-file option
    * Config file specified using the -c my_config option
    * The config file named nextflow.config in the current directory
    * The config file named nextflow.config in the workflow project directory ($projectDir: the directory where the script to be run is located)
    * The config file $HOME/.nextflow/config
    * Values defined within the workflow script itself (e.g., main.nf)

# Resources
[Nextflow's documentation](https://docs.seqera.io/nextflow/)

[Alliancecan's Nextflow wiki](https://docs.alliancecan.ca/wiki/Nextflow/en)

[nf-core pipeline](https://nf-co.re/pipelines)

[Nextflow pipeline example](https://github.com/nextflow-io/awesome-nextflow)

[Nextflow patterns](https://nextflow-io.github.io/patterns/index.html)

[Nextflow containers](https://seqera.io/containers)