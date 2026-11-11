nextflow.enable.dsl = 2

include { CORE } from './workflows/core'
include { DIFFERENTIAL_EXPRESSION } from './workflows/differential_expression'
include { CELL_TYPE_ANALYSIS } from './workflows/cell_type_analysis'
include { INTEGRATION_BENCHMARK } from './workflows/integration_benchmark'
include { TOP_GENE_PROFILES } from './workflows/gene_profiles'
include { AGE_CLOCKS } from './workflows/age_clocks'
include { GENERATE_SYNTHETIC_PSEUDOBULK } from './modules/synthetic_inputs'
include { parseCellTypeSelection } from './modules/cell_type_selection'

workflow {
    def root = file(params.project_dir)
    def pipeline_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def settings = new groovy.json.JsonSlurper().parse(pipeline_file)
    def config = channel.value(pipeline_file)
    def single_cell = null
    def pseudobulk = null
    if (params.analysis_only) {
        single_cell = channel.value(file(params.merged_input, checkIfExists: true))
        pseudobulk = channel.value(file(params.pseudobulk_input, checkIfExists: true))
    } else {
        CORE(root)
        single_cell = CORE.out.single_cell
        pseudobulk = CORE.out.pseudobulk
    }
    def analysis_input = pseudobulk
    if (params.synthetic_analysis) {
        GENERATE_SYNTHETIC_PSEUDOBULK(config)
        analysis_input = GENERATE_SYNTHETIC_PSEUDOBULK.out.pseudobulk
    }
    DIFFERENTIAL_EXPRESSION(
        analysis_input, config,
        channel.value(file(root.resolve('reports/differential_expression_report.py'), checkIfExists: true)),
        channel.value(file(root.resolve('reports/trajectory_analysis_report.py'), checkIfExists: true)),
        channel.value(file(root.resolve('input_data/go'), checkIfExists: true)),
        false, params.synthetic_analysis as boolean, params.allow_partial_results as boolean,
    )
    TOP_GENE_PROFILES(
        analysis_input, config, DIFFERENTIAL_EXPRESSION.out.models,
        params.top_age_gene_count as int, root,
    )
    CELL_TYPE_ANALYSIS(
        single_cell, config,
        channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true)),
        channel.value([
            file(root.resolve(settings.models.aifi_l1 as String), checkIfExists: true),
            file(root.resolve(settings.models.aifi_l2 as String), checkIfExists: true),
        ]),
        parseCellTypeSelection(params.cell_types), DIFFERENTIAL_EXPRESSION.out.models,
    )
    if (params.run_clocks) {
        AGE_CLOCKS(analysis_input, config, params.clock_cell_type ?: '', root)
    }
    if (params.run_benchmark) {
        INTEGRATION_BENCHMARK(
            single_cell, config,
            channel.value(file(root.resolve(settings.models.aifi_l2 as String), checkIfExists: true)),
            settings.integration_benchmark.methods,
            channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true)),
        )
    }
}
