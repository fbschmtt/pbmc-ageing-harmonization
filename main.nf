nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION; clockOptions; profileCountError } from './modules/configuration'

include { CORE } from './workflows/core'
include { DIFFERENTIAL_EXPRESSION } from './workflows/differential_expression'
include { CELL_TYPE_ANALYSIS } from './workflows/cell_type_analysis'
include { INTEGRATION_BENCHMARK } from './workflows/integration_benchmark'
include { TOP_GENE_PROFILES } from './workflows/gene_profiles'
include { AGE_CLOCKS } from './workflows/age_clocks'
include { GENERATE_SYNTHETIC_PSEUDOBULK } from './modules/synthetic_inputs'
include { parseCellTypeSelection } from './modules/cell_type_selection'

workflow {
    def mode = params.mode.toString()
    if (!(mode in ['harmonization', 'analysis', 'full'])) {
        error "Unknown --mode '${mode}'. Choose harmonization, analysis, or full."
    }
    if (mode == 'analysis' && (!params.merged_input || !params.pseudobulk_input)) {
        error "Analysis mode requires both --merged_input and --pseudobulk_input."
    }
    if (mode != 'harmonization') {
        def problem = profileCountError(params.top_age_gene_count)
        if (problem) error(problem)
    }
    def switches = [:]
    [run_benchmark: mode == 'full', run_clocks: false, merge_single_cell: false,
     synthetic_analysis: false, allow_partial_results: false].each { name, default_value ->
        def value = params[name]
        if (value != null && !(value.toString() in ['true', 'false'])) {
            error "--${name} must be true or false."
        }
        // CLI values such as `--run_clocks false` can arrive as strings.
        switches[name] = value == null ? default_value : value.toString().toBoolean()
    }
    def root = file(params.project_dir)
    def pipeline_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def settings = new groovy.json.JsonSlurper().parse(pipeline_file)
    def sections = mode == 'harmonization' ? ['processing'] : ['processing', 'cell_type_analysis', 'differential_expression']
    if (mode != 'harmonization' && switches.run_benchmark) sections << 'integration_benchmark'
    if (mode != 'harmonization' && switches.run_clocks) sections << 'age_clock'
    def registries = mode == 'analysis' ? [] : [settings.studies_config, settings.obs_schema,
        settings.input_sources, 'config/studies.schema.json', 'config/input_sources.schema.json']
        .collect { name -> file(root.resolve(name as String), checkIfExists: true) }
    VALIDATE_CONFIGURATION(channel.value(pipeline_file), sections,
        mode != 'harmonization' && switches.run_clocks ? clockOptions(params) : [:], registries)
    def config = VALIDATE_CONFIGURATION.out.config
    def single_cell = null
    def pseudobulk = null
    if (mode == 'analysis') {
        single_cell = channel.value(file(params.merged_input, checkIfExists: true))
        pseudobulk = channel.value(file(params.pseudobulk_input, checkIfExists: true))
    } else {
        CORE(root, mode == 'full' || switches.merge_single_cell, config)
        single_cell = CORE.out.single_cell
        pseudobulk = CORE.out.pseudobulk
    }
    if (mode != 'harmonization') {
        def analysis_input = pseudobulk
        if (switches.synthetic_analysis) {
            GENERATE_SYNTHETIC_PSEUDOBULK(config)
            analysis_input = GENERATE_SYNTHETIC_PSEUDOBULK.out.pseudobulk
        }
        DIFFERENTIAL_EXPRESSION(
            analysis_input, config,
            channel.value(file(root.resolve('reports/differential_expression_report.py'), checkIfExists: true)),
            channel.value(file(root.resolve('reports/trajectory_analysis_report.py'), checkIfExists: true)),
            channel.value(file(root.resolve('input_data/go'), checkIfExists: true)),
            false, switches.synthetic_analysis, switches.allow_partial_results,
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
        if (switches.run_clocks) {
            AGE_CLOCKS(analysis_input, config, params.clock_cell_type ?: '', root)
        }
        if (switches.run_benchmark) {
            INTEGRATION_BENCHMARK(
                single_cell, config,
                channel.value(file(root.resolve(settings.models.aifi_l2 as String), checkIfExists: true)),
                settings.integration_benchmark.methods,
                channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true)),
            )
        }
    }
}
