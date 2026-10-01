process RENDER_STUDY_QC {
    tag { study_id }
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    tuple val(study_id), path(expression), path(run_report)
    path qc_template

    output:
    path("${study_id}")

    script:
    """
    mkdir -p ${study_id}
    pbmc-qc --project-root ${params.project_dir} --input ${expression} --run-report ${run_report} \\
      --output-dir ${study_id} --study ${study_id} --template ${qc_template}
    """
}

process RENDER_MERGE_QC {
    tag 'all merge outputs'
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    path merged_inputs
    path merge_reports
    path qc_template
    path studies_config

    output:
    path('merged')

    script:
    """
    mkdir -p merged
    pbmc-merge-qc --project-root ${params.project_dir} --input ${merged_inputs} \\
      --run-report ${merge_reports} --output-dir merged --template ${qc_template} \\
      --studies-config ${studies_config}
    """
}
