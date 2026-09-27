process PREPARE_CELLS {
    tag { study_id }
    publishDir "${params.outdir}/prepared", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.cells.csv.gz') || name.endsWith('.source_obs.csv.gz') || name.endsWith('.prepare.json') ? name : null }

    input:
    tuple val(study_id), path(expression)
    path tracked_dependencies, stageAs: 'dependencies/*'

    output:
    tuple val(study_id), path(expression), path("${study_id}.cells.csv.gz"), emit: prepared
    path "${study_id}.source_obs.csv.gz", emit: source_obs

    script:
    """
    pbmc-prepare --project-root ${params.project_dir} --config config/pipeline.json \\
      --study ${study_id} --input ${expression} --output ${study_id}.cells.csv.gz \\
      --source-obs-output ${study_id}.source_obs.csv.gz --report-output ${study_id}.prepare.json
    """
}
