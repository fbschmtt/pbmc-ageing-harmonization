nextflow.enable.dsl = 2

include {
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION
    ANALYSE_CELL_TYPE_DE_DIAGNOSTICS
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS
} from '../../modules/differential_expression'

workflow {
    def pipeline_config = channel.value(file("${launchDir}/config/pipeline.json", checkIfExists: true))
    def fixture = channel.value(file("${projectDir}/fixtures/pseudobulk.txt", checkIfExists: true))
    def go_fixture = channel.value(file("${projectDir}/fixtures", checkIfExists: true))
    def template = channel.value(file("${launchDir}/reports/differential_expression_report.py", checkIfExists: true))
    RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION(
        fixture, pipeline_config, channel.of(tuple('test-type', 'Test type')), false, true,
    )
    ANALYSE_CELL_TYPE_DE_DIAGNOSTICS(
        RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results, pipeline_config, go_fixture,
    )
    RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS(
        RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION.out.results
            .join(ANALYSE_CELL_TYPE_DE_DIAGNOSTICS.out.diagnostics),
        pipeline_config, template,
    )
}
