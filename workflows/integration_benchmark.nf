nextflow.enable.dsl = 2

include { RUN_INTEGRATION_METHOD; AGGREGATE_INTEGRATION_METHODS; RENDER_INTEGRATION_BENCHMARK_REPORT } from '../modules/integration_benchmark'

workflow INTEGRATION_BENCHMARK {
    take:
    merged_input
    pipeline_config
    aifi_l2_model
    method_settings
    report_template

    main:
    def methods = merged_input.flatMap { input ->
        method_settings.collect { method -> tuple(method.name as String, input) }
    }.combine(pipeline_config).combine(aifi_l2_model)
    RUN_INTEGRATION_METHOD(methods)
    AGGREGATE_INTEGRATION_METHODS(
        RUN_INTEGRATION_METHOD.out.status.toList(),
        RUN_INTEGRATION_METHOD.out.artifact.toList(),
        merged_input,
        pipeline_config,
    )
    RENDER_INTEGRATION_BENCHMARK_REPORT(
        AGGREGATE_INTEGRATION_METHODS.out.artifact.combine(AGGREGATE_INTEGRATION_METHODS.out.report),
        report_template,
    )

    emit:
    report = RENDER_INTEGRATION_BENCHMARK_REPORT.out.report
    summary = AGGREGATE_INTEGRATION_METHODS.out.report
}
