process PSEUDOBULK {
    tag { study_id }
    publishDir "${params.outdir}/pseudobulk", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.h5ad') ? name : null }
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    tuple val(study_id), path(expression), path(run_report)
    path merge_config

    output:
    tuple val(study_id), path("${study_id}.pseudobulk.h5ad"), path("${study_id}.pseudobulk.json"), emit: pseudobulk

    script:
    """
    pbmc-merge --project-root . --config ${merge_config} --mode pseudobulk --input ${expression} \\
      --output ${study_id}.pseudobulk.h5ad --report-output ${study_id}.pseudobulk.json
    """
}

process MERGE_PSEUDOBULKS {
    tag 'all studies pseudobulk'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.h5ad') ? name : null }
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    path pseudobulks
    path merge_config

    output:
    tuple val('pseudobulk_merged'), path('pseudobulk_merged.h5ad'), path('pseudobulk_merged.json'), emit: merged

    script:
    """
    pbmc-merge --project-root . --config ${merge_config} --mode pseudobulk-merge --input ${pseudobulks} \\
      --output pseudobulk_merged.h5ad --report-output pseudobulk_merged.json
    """
}

process MERGE_SINGLE_CELLS {
    tag 'all studies single cell'
    publishDir "${params.outdir}/merged", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.h5ad') ? name : null }
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    path studies
    path merge_config
    path aifi_l2_model

    output:
    tuple val('single_cell_merged'), path('single_cell_merged.h5ad'), path('single_cell_merged.json'), emit: merged

    script:
    """
    ${task.ext.harmony_thread_environment ?: ''}
    pbmc-merge --project-root . --config ${merge_config} --mode single-cell-merge \\
      --aifi-l2-model ${aifi_l2_model} --input ${studies} --output single_cell_merged.h5ad \\
      --report-output single_cell_merged.json
    """
}
