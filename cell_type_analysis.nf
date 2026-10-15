nextflow.enable.dsl = 2

include {
    SPLIT_CELL_TYPES
    ANALYSE_CELL_TYPE
    RENDER_CELL_TYPE_REPORT
    WRITE_CELL_TYPE_MANIFEST
    SPLIT_AND_ANALYSE_CELL_TYPES
    ANALYSE_EXISTING_CELL_TYPE_SPLITS
    RENDER_CELL_TYPE_REPORTS
} from './modules/cell_type_analysis'

workflow {
    def root = file(params.project_dir)
    def pipeline_config_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def pipeline_config = channel.value(pipeline_config_file)
    def report_template = channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true))

    def report_inputs = null
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
        report_inputs = primitives.join(analyses)
    } else if (params.cell_type_input) {
        if (params.split_only) {
            error '--cell_type_input and --split_only cannot be used together'
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.parent.name, primitive) }
        ANALYSE_EXISTING_CELL_TYPE_SPLITS(primitives, pipeline_config)
        report_inputs = ANALYSE_EXISTING_CELL_TYPE_SPLITS.out
    } else {
        if (!params.merged_input) {
            error 'Provide --merged_input, for example output/harmonization/merged/single_cell_merged.h5ad'
        }
        def pipeline = new groovy.json.JsonSlurper().parse(pipeline_config_file)
        merged_input = channel.value(file(params.merged_input, checkIfExists: true))
        if (params.split_only) {
            SPLIT_CELL_TYPES(merged_input, pipeline_config)
        } else {
            models = channel.value([
                file(root.resolve(pipeline.models.aifi_l1 as String), checkIfExists: true),
                file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true),
            ])
            SPLIT_AND_ANALYSE_CELL_TYPES(merged_input, pipeline_config)
            report_inputs = SPLIT_AND_ANALYSE_CELL_TYPES.out
        }
    }

    if (report_inputs != null) {
        RENDER_CELL_TYPE_REPORTS(
            report_inputs,
            pipeline_config,
            report_template,
        )
        if (models != null) {
            WRITE_CELL_TYPE_MANIFEST(
                RENDER_CELL_TYPE_REPORTS.out.collect(),
                merged_input,
                pipeline_config,
                models,
            )
        }
    }
}
