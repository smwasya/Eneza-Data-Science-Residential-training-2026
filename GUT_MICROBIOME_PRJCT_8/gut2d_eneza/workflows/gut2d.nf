include { INPUT_CHECK }      from "../modules/local/input_check"
include { FASTQC }           from "../modules/local/quality_control"
include { MULTIQC }          from "../modules/local/quality_control"
include { FILTER_TRIM }      from "../modules/local/filter_trim"
include { DADA2 }            from "../modules/local/dada2"
include { PHYLOSEQ_EXPORT }  from "../modules/local/phyloseq_export"

workflow GUT2D {

    main:

    /*
     * Inputs
     */
    metadata_ch = Channel.fromPath(params.metadata)

    reads_ch = Channel.fromPath(params.reads)

    /*
     * Validate samples
     */
    INPUT_CHECK(
        metadata_ch,
        reads_ch
    )

    /*
     * Read validated sample sheet
     */
    samples_ch = INPUT_CHECK.out.samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row ->

            tuple(
                row.SampleID,
                file(row.forward),
                file(row.reverse)
            )

        }

    /*
     * Duplicate channel
     */
    FASTQC(samples_ch)

MULTIQC(
    FASTQC.out.zip.collect()
)

FILTER_TRIM(samples_ch)
    /*
     * DADA2
     */
    DADA2(

        INPUT_CHECK.out.samplesheet

    )

    /*
     * Phyloseq + ML dataset
     */
    PHYLOSEQ_EXPORT(

        DADA2.out.seqtab,

        DADA2.out.taxonomy,

        metadata_ch

    )

}
