nextflow.enable.dsl = 2

include { GUT2D } from './workflows/gut2d'

workflow {

    GUT2D()

}
