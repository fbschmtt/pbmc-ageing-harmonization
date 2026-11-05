nextflow.enable.dsl = 2

include { RUN_INTEGRATION_METHOD; AGGREGATE_INTEGRATION_METHODS; RENDER_INTEGRATION_BENCHMARK_REPORT } from './modules/integration_benchmark'

workflow {
    def root = file(params.project_dir)
    if (!params.merged_input) {
        error 'Provide --merged_input, for example output/harmonization/merged/single_cell_merged.h5ad'
    }
    def merged_input = file(params.merged_input, checkIfExists: true)
    def pipeline_config = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def document = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
    def aifi_l2_model = file(root.resolve(document.models.aifi_l2 as String), checkIfExists: true)
    def methods = channel.fromList(document.integration_benchmark.methods).map { method ->
        tuple(method.name as String, merged_input, pipeline_config, aifi_l2_model)
    }
    RUN_INTEGRATION_METHOD(methods)
    AGGREGATE_INTEGRATION_METHODS(
        RUN_INTEGRATION_METHOD.out.status.collect(),
        RUN_INTEGRATION_METHOD.out.artifact.filter { artifact -> artifact != null }.collect(),
        channel.value(merged_input),
        channel.value(pipeline_config),
    )
    def report_template = channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true))
    def render_inputs = AGGREGATE_INTEGRATION_METHODS.out.artifact
        .combine(AGGREGATE_INTEGRATION_METHODS.out.report)
    RENDER_INTEGRATION_BENCHMARK_REPORT(render_inputs, report_template)
}
