nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION } from './modules/configuration'

include { CELL_TYPE_ANALYSIS } from './workflows/cell_type_analysis'

include { parseCellTypeSelection; selectCellTypes; pairCellTypeInputs } from './modules/cell_type_selection'

include {
    SPLIT_CELL_TYPES
    ANALYSE_EXISTING_CELL_TYPE_SPLITS
    RENDER_CELL_TYPE_REPORTS
} from './modules/cell_type_analysis'

workflow {
    def selected_slugs = parseCellTypeSelection(params.cell_types)
    if (params.split_only && selected_slugs) {
        error '--cell_types and --split_only cannot be used together'
    }
    def root = file(params.project_dir)
    def pipeline_config_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)

    def report_template = channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true))

    def report_inputs = null
    def analysis_inputs = null
    def merged_input = null
    def models = null
    if (params.cell_type_analysis_input) {
        if (!params.cell_type_input) {
            error '--cell_type_analysis_input requires --cell_type_input'
        }
        if (params.split_only) {
            error '--cell_type_analysis_input and --split_only cannot be used together'
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.parent.name, primitive) }
        def analyses = channel.fromPath(
            params.cell_type_analysis_input, type: 'dir', checkIfExists: true,
        ).map { analysis_dir -> tuple(analysis_dir.parent.name, analysis_dir) }
        report_inputs = primitives.map { slug, path -> tuple('split', slug, path) }
            .mix(analyses.map { slug, path -> tuple('analysis', slug, path) })
            .collect(flat: false)
            .flatMap { entries ->
                def selection = pairCellTypeInputs(entries, selected_slugs)
                if (selection.error) error(selection.error)
                selection.entries
            }
    } else if (params.cell_type_input) {
        if (params.split_only) {
            error '--cell_type_input and --split_only cannot be used together'
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.parent.name, primitive) }
            .collect(flat: false)
            .flatMap { entries ->
                def selection = selectCellTypes(entries, selected_slugs)
                if (selection.error) error(selection.error)
                selection.entries
            }
        analysis_inputs = primitives
    } else {
        if (!params.merged_input) {
            error 'Provide --merged_input, for example output/harmonization/merged/single_cell_merged.h5ad'
        }
        def pipeline = new groovy.json.JsonSlurper().parse(pipeline_config_file)
        merged_input = channel.value(file(params.merged_input, checkIfExists: true))
        if (!params.split_only) models = channel.value([
            file(root.resolve(pipeline.models.aifi_l1 as String), checkIfExists: true),
            file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true),
        ])
    }
    def selected_inputs = report_inputs != null ? report_inputs : analysis_inputs
    def config_input = selected_inputs != null
        ? selected_inputs.toList().map { _inputs -> pipeline_config_file }
        : channel.value(pipeline_config_file)
    VALIDATE_CONFIGURATION(config_input, ['processing', 'cell_type_analysis'], [:], [])
    def pipeline_config = VALIDATE_CONFIGURATION.out.config
    if (analysis_inputs != null) {
        ANALYSE_EXISTING_CELL_TYPE_SPLITS(analysis_inputs, pipeline_config)
        report_inputs = ANALYSE_EXISTING_CELL_TYPE_SPLITS.out
    }
    if (merged_input != null) {
        if (params.split_only) {
            SPLIT_CELL_TYPES(merged_input, pipeline_config)
        } else {
            def de_results = channel.fromPath(
                "${params.outdir}/analysis/per_type/*/differential_expression", type: 'dir', checkIfExists: false,
            ).map { directory -> tuple(directory.parent.name, directory) }
            CELL_TYPE_ANALYSIS(merged_input, pipeline_config, report_template, models, selected_slugs, de_results)
        }
    }

    if (report_inputs != null) {
        def report_inputs_with_de = report_inputs.map { slug, primitive, analysis_dir ->
            def merged_de = file(
                "${params.outdir}/analysis/per_type/${slug}/differential_expression/merged.csv",
                checkIfExists: false,
            )
            def has_merged_de = merged_de.exists()
            def dependency = has_merged_de ? merged_de : pipeline_config_file
            def analysis_files = analysis_dir instanceof List ? analysis_dir
                : (analysis_dir.isDirectory() ? files("${analysis_dir}/*", checkIfExists: true) : [analysis_dir])
            tuple(slug, primitive, analysis_files, has_merged_de, dependency)
        }
        RENDER_CELL_TYPE_REPORTS(
            report_inputs_with_de,
            pipeline_config,
            report_template,
        )

    }
}
