nextflow.enable.dsl = 2

include { INTEGRATION_BENCHMARK } from './workflows/integration_benchmark'

workflow {
    def root = file(params.project_dir)
    if (!params.merged_input) error 'Provide --merged_input'
    def pipeline_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def pipeline = new groovy.json.JsonSlurper().parse(pipeline_file)
    INTEGRATION_BENCHMARK(
        channel.value(file(params.merged_input, checkIfExists: true)),
        channel.value(pipeline_file),
        channel.value(file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true)),
        pipeline.integration_benchmark.methods,
        channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true)),
    )
}
