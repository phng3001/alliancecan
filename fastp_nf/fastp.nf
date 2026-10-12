/*
 * Pipeline parameters
 */

// Primary input

params.reads = "data/*_S*_L*_R{1,2}_*.fastq.gz"
params.outdir = "results"

// Cluster

params.ncpus = 1

/*
 * Check number of CPUs
 */

process NCPUS_CHECK {

    cpus params.ncpus
    
    output:
    stdout

    script:
    """
    echo "Allocated CPUs: ${task.cpus}"
    """
}

/*
 * FastQC before trimming
 */

process FASTQC_BEFORE {

    tag "FastQC before - ${sample_id}"
    publishDir "${params.outdir}/fastqc_before", mode: 'copy'
    // debug true

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.{html,zip}"), emit: reports

    script:
    """
    module load StdEnv/2023 fastqc/0.12.1
    
    fastqc \
        ${reads.join(' ')} \
        -o . \
        -f fastq \
        -q \
        -t ${task.cpus}
    """
}

/*
 * Trimming with FASTP
 */

process FASTP {

    tag "fastp - ${sample_id}"
    publishDir "${params.outdir}/fastp", mode:'copy'
    // debug true

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}/*_paired.fastq.gz"), emit: paired
    tuple val(sample_id), path("${sample_id}/*_unpaired.fastq.gz"), emit: unpaired
    path("${sample_id}/${sample_id}.{html,json}"), emit: reports

    script:
    """
    module load StdEnv/2023 fastp/1.0.1

    mkdir $sample_id

    fastp \\
    -i ${reads[0]} \\
    -I ${reads[1]} \\
    -o ${sample_id}/${sample_id}_R1_paired.fastq.gz \\
    -O ${sample_id}/${sample_id}_R2_paired.fastq.gz \\
    --unpaired1 ${sample_id}/${sample_id}_R1_unpaired.fastq.gz \\
    --unpaired2 ${sample_id}/${sample_id}_R2_unpaired.fastq.gz \\
    --detect_adapter_for_pe \\
    --correction \\
    -w ${task.cpus} \\
    -h ${sample_id}/${sample_id}.html \\
    -j ${sample_id}/${sample_id}.json
    """
}

/*
 * FastQC after trimming
 */

process FASTQC_AFTER {

    tag "FastQC after - ${sample_id}"
    publishDir "${params.outdir}/fastqc_after", mode: 'copy'
    // debug true

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*_fastqc.{html,zip}"), emit: reports

    script:
    """
    module load StdEnv/2023 fastqc/0.12.1

    fastqc \
        ${reads.join(' ')} \
        -o . \
        -f fastq \
        -q \
        -t ${task.cpus}
    """
}

/*
 * MultiQC
 */

process MULTIQC {
    container "/project/def-mouellet/Scripts_MOU/PNP/containers/multiqc_v1.25.1.sif"
    tag "MultiQC"
    publishDir "${params.outdir}/multiqc", mode: 'copy'
    // debug true

    input:
    path('*')

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc . --interactive
    """
}



workflow {
    // Workflow parameters
    log.info"""\
        FASTP NF PIPELINE
        ===================================
        reads : ${params.reads}
        outdir: ${params.outdir}
        ===================================
        """
        .stripIndent()

    // Input chanel
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists:true)
    read_pairs_ch.view()

    // Check number of CPUs
    NCPUS_CHECK().view()

    // Raw read FastQC
    FASTQC_BEFORE(read_pairs_ch)

    // Trimming with FASTP
    FASTP(read_pairs_ch)

    // Trimmed read FastQC
    FASTQC_AFTER(FASTP.out.paired)

    // MultiQC
    multiqc_ch = FASTQC_BEFORE.out.reports
        .mix(FASTQC_AFTER.out.reports)
        .mix(FASTP.out.reports)
        .collect()
    // multiqc_ch.view()
    MULTIQC(multiqc_ch)

    // Workflow completion handler
    workflow.onComplete = {
    if (workflow.success) {
        log.info"""
        ===================================
        Pipeline completed successfully!
        MultiQC report:
        ${params.outdir}/multiqc/multiqc_report.html
        ===================================
        """.stripIndent()
    } else {
        log.error("Pipeline failed. Check the Nextflow log for details.")
        }
    }
}
