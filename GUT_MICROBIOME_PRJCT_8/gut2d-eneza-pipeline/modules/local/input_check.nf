process INPUT_CHECK {

    tag "Input validation"

    publishDir "${params.outdir}/validation", mode: "copy"

    //container "bioconductor/bioconductor_docker:RELEASE_3_22"

    input:
    path metadata
    path reads

    output:
    path "validated_samplesheet.tsv", emit: samplesheet

    script:
    """
    Rscript ${projectDir}/bin/input_check.R \
        --metadata ${metadata} \
        --reads ${reads} \
        --output validated_samplesheet.tsv
    """
}
