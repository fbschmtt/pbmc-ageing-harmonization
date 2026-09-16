process WRITE_RUN_MANIFEST {
    tag 'published artifacts'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path artifacts
    path reports
    path pipeline_config
    val selected_studies

    output:
    path 'run_manifest.json'

    script:
    """
    pbmc-manifest --config ${pipeline_config} --artifact ${artifacts} --report ${reports} \\
      --studies '${selected_studies}' --pipeline-revision '${params.pipeline_revision}' \\
      --python-image '${params.python_image}' --r-image '${params.r_image}' --output run_manifest.json
    """
}
