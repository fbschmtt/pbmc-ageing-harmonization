process RUN_INTEGRATION_BENCHMARK {
    tag 'global integration benchmark'
    publishDir "${params.outdir}/integration_benchmark", mode: 'copy', overwrite: true

    input:
    path merged_input
    path pipeline_config
    path aifi_l2_model

    output:
    path 'integration_benchmark.h5ad', emit: artifact
    path 'integration_benchmark.json', emit: report

    script:
    """
    pbmc-integration-benchmark --project-root . --config ${pipeline_config} \\
      --input ${merged_input} --aifi-l2-model ${aifi_l2_model} \\
      --output integration_benchmark.h5ad --report-output integration_benchmark.json
    """
}

process RENDER_INTEGRATION_BENCHMARK_REPORT {
    tag 'global integration benchmark report'
    publishDir "${params.outdir}/integration_benchmark", mode: 'copy', overwrite: true

    input:
    path artifact
    path run_report
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
