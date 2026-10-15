process LIST_PSEUDOBULK_CELL_TYPES {
    tag 'list pseudobulk cell types'

    input:
    path pseudobulk_input
    path pipeline_config

    output:
    path 'cell_types.tsv', emit: cell_types

    script:
    """
    pbmc-differential-expression --input ${pseudobulk_input} --config ${pipeline_config} \
      --list-cell-types --cell-types-output cell_types.tsv
    """
}

process RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION {
    tag { cell_type }
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            filename == "cell_type_result/${slug}"
                ? "${slug}/differential_expression"
                : null
        }

    input:
    path pseudobulk_input
    path pipeline_config
    tuple val(slug), val(cell_type)
    val de_test_mode
    val synthetic_test_data

    output:
    path 'cell_type_result/*', emit: results

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir cell_type_result --config ${pipeline_config} --cell-type "${cell_type}" \
      --cpus ${task.cpus} ${test_mode_argument} ${synthetic_data_argument}
    """
}

process WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST {
    tag 'collect cell-type PyDESeq2 models'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'analysis/per_type' ? null : filename }

    input:
    path pseudobulk_input
    path pipeline_config
    path cell_type_results
    val de_test_mode
    val synthetic_test_data

    output:
    path 'analysis/combined/differential_expression', emit: results
    path 'analysis/per_type', emit: per_type

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir analysis/combined/differential_expression \
      --per-type-output-dir analysis/per_type --config ${pipeline_config} \
      --combine-cell-type-results ${cell_type_results} \
      ${test_mode_argument} ${synthetic_data_argument}
    """
}

process RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT {
    tag 'cross-cell-type trajectory and expression analysis'
    publishDir "${params.outdir}/analysis/combined", mode: 'copy', overwrite: true

    input:
    path per_type_de
    path pseudobulk_input
    path pipeline_config
    path report_template

    output:
    path 'trajectory', emit: report

    script:
    """
    pbmc-trajectory-report --differential-expression-dir ${per_type_de} \\
      --pseudobulk-input ${pseudobulk_input} --output-dir trajectory --config ${pipeline_config} \\
      --template ${report_template} --project-root .
    """
}

process RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS {
    tag 'per-cell-type differential expression reports'
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            filename.startsWith('de_reports/') ? filename.substring('de_reports/'.length()) : null
        }

    input:
    path per_type_de
    path cross_cell_type_report

    output:
    path 'de_reports/*/differential_expression/report.html', emit: reports

    script:
    """
    test -s ${cross_cell_type_report}/report.html
    pbmc-differential-expression-report --input-dir ${per_type_de} --output-root de_reports
    """
}
