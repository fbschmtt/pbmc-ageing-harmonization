nextflow.enable.dsl = 2

include {
    LIST_PSEUDOBULK_CELL_TYPES
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION
    WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST
    RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS
    ANALYSE_CELL_TYPE_DE_DIAGNOSTICS
} from '../modules/differential_expression'

workflow DIFFERENTIAL_EXPRESSION {
    take:
    pseudobulk_input
    pipeline_config
    de_report_template
    trajectory_report_template
    go_enrichment_data
    de_test_mode
    synthetic_test_data
    allow_partial_results

    main:
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
            .toList(),
        de_test_mode,
        synthetic_test_data,
        allow_partial_results,
    )
    RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT(
        WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST.out.per_type,
        pseudobulk_input,
        pipeline_config,
        trajectory_report_template,
    )
    ANALYSE_CELL_TYPE_DE_DIAGNOSTICS(
        RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results,
        pipeline_config,
        go_enrichment_data,
    )
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS(
        RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results
            .join(ANALYSE_CELL_TYPE_DE_DIAGNOSTICS.out.diagnostics),
        pipeline_config,
        de_report_template,
    )

    emit:
    models = RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results
    per_type = WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST.out.per_type
    combined = WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST.out.results
    reports = RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS.out.reports
    cross_type = RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT.out.report
}
