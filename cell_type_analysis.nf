nextflow.enable.dsl = 2

include { SPLIT_CELL_TYPES; ANALYSE_CELL_TYPE; RENDER_CELL_TYPE_REPORT; WRITE_CELL_TYPE_MANIFEST } from './modules/cell_type_analysis'

workflow {
    def root = file(params.project_dir)
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def report_template = channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true))

    if (params.cell_type_analysis_input) {
        if (!params.cell_type_input) {
            error "--cell_type_analysis_input requires --cell_type_input"
        }
        if (params.split_only) {
            error "--cell_type_analysis_input and --split_only cannot be used together"
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.baseName, primitive) }
        def analyses = channel.fromPath(
            params.cell_type_analysis_input, type: 'dir', checkIfExists: true,
        ).map { analysis_dir -> tuple(analysis_dir.baseName, analysis_dir) }
        RENDER_CELL_TYPE_REPORT(primitives.join(analyses), pipeline_config, report_template)
    } else if (params.cell_type_input) {
        if (params.split_only) {
            error "--cell_type_input and --split_only cannot be used together"
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.baseName, primitive) }
        ANALYSE_CELL_TYPE(primitives, pipeline_config)
        RENDER_CELL_TYPE_REPORT(
            primitives.join(ANALYSE_CELL_TYPE.out.analysis), pipeline_config, report_template,
        )
    } else {
        if (!params.merged_input) {
            error "Provide --merged_input, for example output/merged/single_cell_merged.h5ad"
        }
        def pipeline = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
        def merged_input = channel.value(file(params.merged_input, checkIfExists: true))
        SPLIT_CELL_TYPES(merged_input, pipeline_config)
        if (!params.split_only) {
            def models = channel.value([
                file(root.resolve(pipeline.models.aifi_l1 as String), checkIfExists: true),
                file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true),
            ])
            def primitives = SPLIT_CELL_TYPES.out.cells.flatten()
                .map { primitive -> tuple(primitive.baseName, primitive) }
            ANALYSE_CELL_TYPE(primitives, pipeline_config)
            RENDER_CELL_TYPE_REPORT(
                primitives.join(ANALYSE_CELL_TYPE.out.analysis), pipeline_config, report_template,
            )
            WRITE_CELL_TYPE_MANIFEST(RENDER_CELL_TYPE_REPORT.out.collect(), merged_input, pipeline_config, models)
        }
    }
}
