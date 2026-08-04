process FASTQC {

    tag "${sample}"

    publishDir "${params.outdir}/qc/fastqc", mode: "copy"

   // container "biocontainers/fastqc:v0.12.1_cv8"

    cpus 4

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    path "*_fastqc.zip", emit: zip
    path "*_fastqc.html", emit: html

    script:
    """
    fastqc \
        --threads ${task.cpus} \
        ${read1} \
        ${read2}
    """
}

process MULTIQC {

    publishDir "${params.outdir}/qc", mode: "copy"

   // container "ewels/multiqc:1.32"

    input:
    path reports

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data

    script:
    """
    multiqc . --force
    """
}
