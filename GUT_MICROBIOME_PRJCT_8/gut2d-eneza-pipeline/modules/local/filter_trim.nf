process FILTER_TRIM {

    tag "${sample}"

    publishDir "${params.outdir}/fastp", mode: "copy"

    cpus 3

    input:
    tuple val(sample), path(read1), path(read2)

    output:
    tuple val(sample),
          path("${sample}_R1.trimmed.fastq.gz"),
          path("${sample}_R2.trimmed.fastq.gz"),
          emit: trimmed_reads

    path "${sample}.html", emit: html
    path "${sample}.json", emit: json

    script:
    """
    fastp \
        -i ${read1} \
        -I ${read2} \
        -o ${sample}_R1.trimmed.fastq.gz \
        -O ${sample}_R2.trimmed.fastq.gz \
        -h ${sample}.html \
        -j ${sample}.json \
        --thread ${task.cpus}
    """
}
