nextflow.enable.dsl = 2

include { RUN_PSEUDOBULK_DIFFERENTIAL_EXPRESSION } from './modules/differential_expression'

workflow {
    def root = file(params.project_dir)
    if (!params.pseudobulk_input) {
        error 'Provide --pseudobulk_input, for example output/merged/pseudobulk_merged.h5ad'
    }
    def pseudobulk_input = channel.value(file(params.pseudobulk_input, checkIfExists: true))
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def de_test_mode = channel.value(params.de_test_mode as boolean)
    def synthetic_test_data = channel.value(params.de_synthetic_test_data as boolean)
    RUN_PSEUDOBULK_DIFFERENTIAL_EXPRESSION(
        pseudobulk_input,
        pipeline_config,
        de_test_mode,
        synthetic_test_data,
    )
}
