nextflow.enable.dsl = 2

include {
    LIST_PSEUDOBULK_CELL_TYPES
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION
    WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST
    RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS
} from './modules/differential_expression'

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
        def report_inputs = channel.fromPath(
            "${params.per_type_de_dir}/*/differential_expression",
            type: 'dir',
            checkIfExists: true,
        ).map { result_dir -> tuple(result_dir.parent.name, result_dir) }
        RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS(
            report_inputs,
            pipeline_config,
            de_report_template,
            go_enrichment_data,
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
            RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results
                .map { _slug, result_dir -> result_dir }
                .collect(),
            de_test_mode,
            synthetic_test_data,
        )
        RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT(
            WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST.out.per_type,
            pseudobulk_input,
            pipeline_config,
            trajectory_report_template,
        )
        RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS(
            RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results,
            pipeline_config,
            de_report_template,
            go_enrichment_data,
        )
    }
}
