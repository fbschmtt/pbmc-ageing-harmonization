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

    stub:
    // The wildcard fails whichever method the executor starts first, independent of scheduling order.
    def fail = params.stub_fail_benchmark_method == '*' ? """
    if mkdir '${params.outdir}/.benchmark-failed-once' 2>/dev/null; then exit 53; fi
    """ : (method_name == params.stub_fail_benchmark_method ? 'exit 53' : '')
    """
    ${fail}
    printf '%s\n' 'Benchmark stub' > method_${method_name}.h5ad
    printf '%s\n' '{"name":"${method_name}","status":"complete"}' > method_${method_name}.json
    """
}

process AGGREGATE_INTEGRATION_METHODS {
    errorStrategy 'ignore'
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

    stub:
    """
    python3 - <<'PY_STUB'
    import json
    from pathlib import Path
    expected = {method['name'] for method in json.loads(Path('${pipeline_config}').read_text())['integration_benchmark']['methods']}
    methods = [json.loads(path.read_text()) for path in Path('statuses').glob('*.json')]
    missing = sorted(expected - {method['name'] for method in methods})
    Path('integration_benchmark.json').write_text(json.dumps({
        'status': 'partial' if missing else 'complete', 'methods': methods, 'failed_methods': missing,
    }))
    Path('integration_benchmark.h5ad').write_text('Benchmark aggregate stub')
    PY_STUB
    """
}

process RENDER_INTEGRATION_BENCHMARK_REPORT {
    errorStrategy 'ignore'
    tag 'global integration benchmark report'
    publishDir "${params.outdir}/harmonization/integration_benchmark", mode: 'copy', overwrite: true

    input:
    tuple path(artifact), path(run_report)
    path report_template

    output:
    path 'report.html', emit: report
    path 'executed.ipynb', emit: notebook

    script:
    """
    pbmc-integration-benchmark-report --project-root . --input ${artifact} \\
      --run-report ${run_report} --template ${report_template} --output-dir .
    """

    stub:
    """
    test -f ${artifact}
    test -f ${run_report}
    printf '%s\n' '<html>Benchmark report stub</html>' > report.html
    printf '%s\n' '{}' > executed.ipynb
    """
}
