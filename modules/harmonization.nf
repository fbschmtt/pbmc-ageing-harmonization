process HARMONIZE {
    tag { study_id }
    publishDir "${params.outdir}/harmonized", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.h5ad') ? name : null }
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    tuple val(study_id), path(expression), path(prepared_cells)
    path tracked_dependencies, stageAs: 'dependencies/*'

    output:
    tuple val(study_id), path("${study_id}.h5ad"), path("${study_id}.json"), emit: harmonized

    script:
    """
    pbmc-harmonize --project-root ${params.project_dir} --config config/pipeline.json \\
      --study ${study_id} --input ${expression} --prepared-obs ${prepared_cells} \\
      --output ${study_id}.h5ad --report-output ${study_id}.json
    """
}
