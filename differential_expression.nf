nextflow.enable.dsl = 2

include {
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS
    ANALYSE_CELL_TYPE_DE_DIAGNOSTICS
} from './modules/differential_expression'

include { DIFFERENTIAL_EXPRESSION } from './workflows/differential_expression'

workflow {
    def root = file(params.project_dir)
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def de_report_template = channel.value(
        file(root.resolve('reports/differential_expression_report.py'), checkIfExists: true)
    )
    def go_enrichment_data = channel.value(
        file(root.resolve('input_data/go'), checkIfExists: true)
    )
    if (params.per_type_de_dir) {
        def report_glob = params.de_cell_type
            ? "${params.per_type_de_dir}/${params.de_cell_type}/differential_expression"
            : "${params.per_type_de_dir}/*/differential_expression"
        def report_inputs = channel.fromPath(
            report_glob,
            type: 'dir',
            checkIfExists: true,
        ).map { result_dir -> tuple(result_dir.parent.name, result_dir) }
        ANALYSE_CELL_TYPE_DE_DIAGNOSTICS(report_inputs, pipeline_config, go_enrichment_data)
        RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS(
            report_inputs.join(ANALYSE_CELL_TYPE_DE_DIAGNOSTICS.out.diagnostics),
            pipeline_config,
            de_report_template,
        )
    } else {
        if (!params.pseudobulk_input) {
            error 'Provide --pseudobulk_input, for example output/harmonization/merged/pseudobulk_merged.h5ad'
        }
        def pseudobulk_input = channel.value(file(params.pseudobulk_input, checkIfExists: true))
        def trajectory_report_template = channel.value(
            file(root.resolve('reports/trajectory_analysis_report.py'), checkIfExists: true)
        )
        def de_test_mode = channel.value(params.de_test_mode as boolean)
        def synthetic_test_data = channel.value(params.de_synthetic_test_data as boolean)
        def allow_partial_results = channel.value(params.allow_partial_results as boolean)
        DIFFERENTIAL_EXPRESSION(
            pseudobulk_input, pipeline_config, de_report_template, trajectory_report_template,
            go_enrichment_data, de_test_mode, synthetic_test_data, allow_partial_results,
        )
    }
}
