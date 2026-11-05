process RUN_INTEGRATION_METHOD {
    tag { "integration method ${method_name}" }
    maxForks 1
    errorStrategy 'ignore'
    publishDir "${params.outdir}/harmonization/integration_benchmark/methods", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    tuple val(method_name), path(merged_input), path(pipeline_config), path(aifi_l2_model)

    output:
    path "method_${method_name}.h5ad", emit: artifact, optional: true
    path "method_${method_name}.json", emit: status

    script:
    """
    pbmc-integration-benchmark --project-root . --config ${pipeline_config} \\
      --input ${merged_input} --log-level ${params.log_level} method \\
      --method-name ${method_name} --aifi-l2-model ${aifi_l2_model} \\
      --output method_${method_name}.h5ad --status-output method_${method_name}.json
    """
}

process AGGREGATE_INTEGRATION_METHODS {
    tag 'aggregate integration benchmark methods'
    publishDir "${params.outdir}/harmonization/integration_benchmark", mode: 'copy', overwrite: true

    input:
    path method_statuses, stageAs: 'statuses/*'
    path method_artifacts, stageAs: 'artifacts/*'
    path merged_input
    path pipeline_config

    output:
    path 'integration_benchmark.h5ad', emit: artifact, optional: true
    path 'integration_benchmark.json', emit: report

    script:
    """
    pbmc-integration-benchmark --project-root . --config ${pipeline_config} \\
      --input ${merged_input} --log-level ${params.log_level} aggregate \\
      ${method_statuses ? "--method-statuses ${method_statuses.join(' ')}" : ''} \\
      ${method_artifacts ? "--method-artifacts ${method_artifacts.join(' ')}" : ''} \\
      --output integration_benchmark.h5ad --report-output integration_benchmark.json
    """
}

process RENDER_INTEGRATION_BENCHMARK_REPORT {
    tag 'global integration benchmark report'
    publishDir "${params.outdir}/harmonization/integration_benchmark", mode: 'copy', overwrite: true

    input:
    tuple path(artifact), path(run_report)
    path report_template

    output:
    path 'report.html'
    path 'executed.ipynb'

    script:
    """
    pbmc-integration-benchmark-report --project-root . --input ${artifact} \\
      --run-report ${run_report} --template ${report_template} --output-dir .
    """
}
