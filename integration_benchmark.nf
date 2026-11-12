nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION } from './modules/configuration'

include { INTEGRATION_BENCHMARK } from './workflows/integration_benchmark'

workflow {
    def root = file(params.project_dir)
    if (!params.merged_input) error 'Provide --merged_input'
    def pipeline_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def pipeline = new groovy.json.JsonSlurper().parse(pipeline_file)
    VALIDATE_CONFIGURATION(channel.value(pipeline_file), ['processing', 'integration_benchmark'], [:], [])
    INTEGRATION_BENCHMARK(
        channel.value(file(params.merged_input, checkIfExists: true)),
        VALIDATE_CONFIGURATION.out.config,
        channel.value(file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true)),
        pipeline.integration_benchmark.methods,
        channel.value(file(root.resolve('reports/integration_benchmark_report.py'), checkIfExists: true)),
    )
}
