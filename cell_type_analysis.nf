nextflow.enable.dsl = 2

include { SPLIT_CELL_TYPES; RUN_CELL_TYPE_REPORT; WRITE_CELL_TYPE_MANIFEST } from './modules/cell_type_analysis'

workflow {
    if (!params.merged_input) {
        error "Provide --merged_input, for example output/merged/single_cell_merged.h5ad"
    }
    def root = file(params.project_dir)
    def pipeline = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
    def merged_input = channel.value(file(params.merged_input, checkIfExists: true))
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def models = channel.value([
        file(root.resolve(pipeline.models.aifi_l1 as String), checkIfExists: true),
        file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true),
    ])
    def report_template = channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true))
    SPLIT_CELL_TYPES(merged_input, pipeline_config)
    RUN_CELL_TYPE_REPORT(SPLIT_CELL_TYPES.out.cells.flatten(), pipeline_config, report_template)
    WRITE_CELL_TYPE_MANIFEST(RUN_CELL_TYPE_REPORT.out.collect(), merged_input, pipeline_config, models)
}
