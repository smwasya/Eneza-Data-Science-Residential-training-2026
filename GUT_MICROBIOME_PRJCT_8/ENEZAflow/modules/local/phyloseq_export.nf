process PHYLOSEQ_EXPORT {

    tag "Phyloseq"

    publishDir "${params.outdir}/phyloseq",
        mode: "copy"

    input:

    path seqtab

    path taxonomy

    path metadata

    output:

    path "phyloseq.rds"

    path "phyloseq_genus.rds"

    path "machine_learning_dataset.csv"

    path "taxonomy_table.csv"

    path "asv_abundance_table.csv"

    path "rep_seqs.fasta"

    script:

    """
    phyloseq_export.R \
        --seqtab ${seqtab} \
        --taxonomy ${taxonomy} \
        --metadata ${metadata} \
        --output .
    """
}
