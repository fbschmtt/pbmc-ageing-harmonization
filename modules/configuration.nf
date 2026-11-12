def profileCountError(value) {
    def raw = value?.toString()
    raw && raw.isInteger() && raw.toInteger() > 0 ? null : '--top_age_gene_count must be a positive integer'
}

def clockOptions(params) {
    [external_studies: params.clock_external_studies, methods: params.clock_methods,
     alphas: params.clock_alphas, l1_ratios: params.clock_l1_ratios, n_jobs: params.clock_n_jobs,
     calibration_samples: params.clock_calibration_samples, calibration_seed: params.clock_calibration_seed]
}

process VALIDATE_CONFIGURATION {
    tag 'configuration preflight'
    errorStrategy 'terminate'
    container params.python_image

    input:
    path pipeline_config
    val sections
    val clock_options
    path registry_files, stageAs: 'configuration/*'

    output:
    path 'validated/pipeline.json', emit: config
    path 'validated/preflight.json', emit: validation

    script:
    def clock_json = groovy.json.JsonOutput.toJson(clock_options).replace("'", "'\"'\"'")
    def core_arg = registry_files ? "--project-root '${params.project_dir}'" : ''
    """
    python -m pbmc_pipeline.preflight --config ${pipeline_config} --sections '${sections.join(',')}' \\
      --clock-options '${clock_json}' ${core_arg}
    """

    stub:
    def fail = params.stub_fail_preflight ? 'exit 61' : ''
    """
    ${fail}
    mkdir validated
    cp ${pipeline_config} validated/pipeline.json
    printf '%s\n' '{"status":"valid","test_only":true}' > validated/preflight.json
    """
}
