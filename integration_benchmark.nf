nextflow.enable.dsl = 2

include { RUN_INTEGRATION_BENCHMARK; RENDER_INTEGRATION_BENCHMARK_REPORT } from './modules/integration_benchmark'

workflow {
    def root = file(params.project_dir)
    if (!params.merged_input) {
        error 'Provide --merged_input, for example output/merged/single_cell_merged.h5ad'
    }
    def merged_input = channel.value(file(params.merged_input, checkIfExists: true))
    def pipeline_config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def document = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
    def aifi_l2_model = channel.value(file(root.resolve(document.models.aifi_l2 as String), checkIfExists: true))
    RUN_INTEGRATION_BENCHMARK(merged_input, pipeline_config, aifi_l2_model)
    def report_template = channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true))
    RENDER_INTEGRATION_BENCHMARK_REPORT(
        RUN_INTEGRATION_BENCHMARK.out.artifact,
        RUN_INTEGRATION_BENCHMARK.out.report,
        report_template,
    )
}
