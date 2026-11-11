nextflow.enable.dsl = 2

include { CORE } from './workflows/core'

workflow {
    CORE(file(params.project_dir))
}
