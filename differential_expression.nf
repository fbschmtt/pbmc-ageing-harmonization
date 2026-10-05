nextflow.enable.dsl = 2

include {
    LIST_PSEUDOBULK_CELL_TYPES
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION
    WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST
    RENDER_CROSS_CELL_TYPE_TRAJECTORY_REPORT
} from './modules/differential_expression'

workflow {
    def root = file(params.project_dir)
    if (!params.pseudobulk_input) {
        error 'Provide --pseudobulk_input, for example output/merged/pseudobulk_merged.h5ad'
    }
    def pseudobulk_input = channel.value(file(params.pseudobulk_input, checkIfExists: true))
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def trajectory_report_template = channel.value(
        file(root.resolve('reports/trajectory_analysis_report.py'), checkIfExists: true)
    )
    def de_test_mode = channel.value(params.de_test_mode as boolean)
    def synthetic_test_data = channel.value(params.de_synthetic_test_data as boolean)
    LIST_PSEUDOBULK_CELL_TYPES(pseudobulk_input, pipeline_config)
    def cell_types = LIST_PSEUDOBULK_CELL_TYPES.out.cell_types
        .splitText()
        .map { line ->
            def fields = line.trim().split('\\t', 2)
            if (fields.size() != 2 || !fields[0] || !fields[1]) {
                error "Invalid pseudobulk cell-type index row: ${line.inspect()}"
            }
            tuple(fields[0], fields[1])
        }
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION(
        pseudobulk_input,
        pipeline_config,
        cell_types,
        de_test_mode,
        synthetic_test_data,
    )
    WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST(
        pseudobulk_input,
        pipeline_config,
        RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results.collect(),
        de_test_mode,
        synthetic_test_data,
    )
    RENDER_CROSS_CELL_TYPE_TRAJECTORY_REPORT(
        WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST.out.results,
        pipeline_config,
        trajectory_report_template,
    )
}
