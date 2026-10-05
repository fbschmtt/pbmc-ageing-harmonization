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
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path pseudobulk_input
    path pipeline_config
    path cell_type_results
    val de_test_mode
    val synthetic_test_data

    output:
    path 'differential_expression', emit: results

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir differential_expression --config ${pipeline_config} \
      --combine-cell-type-results ${cell_type_results} \
      ${test_mode_argument} ${synthetic_data_argument}
    """
}

process RENDER_CROSS_CELL_TYPE_TRAJECTORY_REPORT {
    tag 'cross-cell-type age trajectory report'
    publishDir "${params.outdir}/differential_expression", mode: 'copy', overwrite: true

    input:
    path differential_expression
    path pipeline_config
    path report_template

    output:
    path 'trajectory_analysis', emit: report

    script:
    """
    pbmc-trajectory-report --differential-expression-dir ${differential_expression} \
      --output-dir trajectory_analysis --config ${pipeline_config} \
      --template ${report_template} --project-root .
    """
}
