process LIST_CLOCK_CELL_TYPES {
    input:
    path pseudobulk_input
    path pipeline_config
    val selection
    path index_script

    output:
    path 'cell_types.tsv', emit: types
    path 'expected.json', emit: expected

    script:
    """
    python ${index_script} --input ${pseudobulk_input} --config ${pipeline_config} \\
      --selection '${selection.replace("'", "'\"'\"'")}' --output cell_types.tsv --expected expected.json
    """

    stub:
    """
    printf 'test-type\tTest type\nfailed-type\tFailed type\n' > cell_types.tsv
    printf '%s\n' '["Test type", "Failed type"]' > expected.json
    """
}

process FIT_AGE_CLOCK {
    tag { slug }
    errorStrategy 'ignore'
    maxForks 1
    publishDir { params.age_clock_output ?: "${params.outdir}/analysis/per_type" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == "${slug}/ageing_clock" ? filename : null }
    publishDir "${params.outdir}/analysis/combined/age_clocks", mode: 'copy', overwrite: true,
        pattern: '*.status.json'

    input:
    path pseudobulk_input
    path pipeline_config
    tuple val(slug), val(cell_type)
    val input_identity
    path clock_script

    output:
    tuple val(slug), path("${slug}/ageing_clock"), emit: model, optional: true
    path "${slug}.status.json", emit: status

    script:
    """
    python ${clock_script} fit --input ${pseudobulk_input} --config ${pipeline_config} --output-dir . \\
      --cell-type '${cell_type.replace("'", "'\"'\"'")}' --slug ${slug} --status ${slug}.status.json \\
      --input-identity '${input_identity.replace("'", "'\"'\"'")}' \\
      --external-studies '${params.clock_external_studies.toString().replace("'", "'\"'\"'")}' --methods '${(params.clock_methods ?: '').toString().replace("'", "'\"'\"'")}' \\
      --alphas '${params.clock_alphas.toString().replace("'", "'\"'\"'")}' --l1-ratios '${params.clock_l1_ratios.toString().replace("'", "'\"'\"'")}' \\
      --n-jobs ${params.clock_n_jobs as int} --calibration-samples ${params.clock_calibration_samples as int} --calibration-seed ${params.clock_calibration_seed as int}
    """

    stub:
    def fail = params.stub_fail_clock && slug == 'failed-type' ? 'exit 19' : ''
    """
    ${fail}
    mkdir -p ${slug}/ageing_clock/elastic_net
    printf '%s\n' '{"cell_type":"${cell_type}","methods":[{"method":"elastic_net"}]}' > ${slug}/ageing_clock/age_clock.json
    printf '%s\n' '{"cell_type":"${cell_type}","slug":"${slug}","status":"complete"}' > ${slug}.status.json
    """
}

process RENDER_AGE_CLOCK {
    tag { slug }
    errorStrategy 'ignore'
    publishDir { params.age_clock_output ?: "${params.outdir}/analysis/per_type" }, mode: 'copy', overwrite: true

    input:
    tuple val(slug), path(model_dir, stageAs: 'prepared')
    path method_template
    path comparison_template

    output:
    path "${slug}/ageing_clock/{report.html,executed.ipynb,*/report.html,*/executed.ipynb}", emit: reports

    script:
    """
    mkdir -p ${slug}/ageing_clock
    cp -a prepared/. ${slug}/ageing_clock/
    python -m pbmc_pipeline.age_clock_report --report-root . --cell-type ${slug} \\
      --template ${method_template} --comparison-template ${comparison_template} --project-root .
    """

    stub:
    def fail = params.stub_fail_clock_report && slug == 'failed-type' ? 'exit 37' : ''
    """
    ${fail}
    mkdir -p ${slug}/ageing_clock/elastic_net
    printf '%s\n' '<html>Clock comparison stub</html>' > ${slug}/ageing_clock/report.html
    printf '%s\n' '{}' > ${slug}/ageing_clock/executed.ipynb
    printf '%s\n' '<html>Clock method stub</html>' > ${slug}/ageing_clock/elastic_net/report.html
    printf '%s\n' '{}' > ${slug}/ageing_clock/elastic_net/executed.ipynb
    """
}

process AGGREGATE_AGE_CLOCKS {
    errorStrategy 'ignore'
    publishDir { params.age_clock_output ?: "${params.outdir}/analysis/per_type" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.startsWith('types/') ? filename.substring('types/'.length()) : filename }

    input:
    tuple val(slugs), path(model_dirs, stageAs: 'models/model??')
    path statuses, stageAs: 'statuses/*'
    val expected_labels
    val existing
    path clock_script
    path cross_type_template

    output:
    path 'types/age_clock_cross_type', emit: aggregate
    path 'age_clock_all_types_summary.json', emit: summary

    script:
    def expected_json = groovy.json.JsonOutput.toJson(expected_labels)
    def directories = model_dirs instanceof List ? model_dirs : [model_dirs]
    def stage_models = directories.withIndex().collect { directory, index ->
        "mkdir -p types/${slugs[index]}/ageing_clock; cp -a ${directory}/. types/${slugs[index]}/ageing_clock/"
    }.join('\n')
    """
    mkdir -p types statuses
    ${stage_models}
    printf '%s\n' '${expected_json.replace("'", "'\"'\"'")}' > expected.json
    python ${clock_script} aggregate --input-root types --status-dir statuses --expected expected.json \\
      --output age_clock_all_types_summary.json --template ${cross_type_template} ${existing ? '--existing' : ''}
    """

    stub:
    """
    mkdir -p types/age_clock_cross_type
    printf '%s\n' '{"methods":[]}' > types/age_clock_cross_type/summary.json
    printf '%s\n' '{"test_only":true}' > age_clock_all_types_summary.json
    """
}
