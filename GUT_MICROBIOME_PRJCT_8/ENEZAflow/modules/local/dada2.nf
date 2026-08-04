process DADA2 {

    tag "DADA2"

    publishDir "${params.outdir}/dada2",
        mode: "copy"

    cpus params.dada2.cpus

    memory params.dada2.memory

    input:

    path samplesheet

    output:

    path "results/seqtab.rds", emit: seqtab_raw

    path "results/seqtab_nochim.rds", emit: seqtab

    path "results/taxonomy.rds", emit: taxonomy

    path "results/taxonomy.csv"

    path "results/error_rates.pdf"

    path "results/filtering_summary.csv"

    script:

    """
    dada2_pipeline.R \
        --samplesheet ${samplesheet} \
        --silva ${params.silva} \
        --threads ${task.cpus} \
        --output results
    """
}
