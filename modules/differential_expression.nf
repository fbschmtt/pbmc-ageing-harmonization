process RUN_PSEUDOBULK_DIFFERENTIAL_EXPRESSION {
    tag 'per-study and merged PyDESeq2 models'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path pseudobulk_input
    path pipeline_config
    val de_test_mode
    val synthetic_test_data

    output:
    path 'differential_expression', emit: results

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir differential_expression --config ${pipeline_config} --cpus ${task.cpus} \
      ${test_mode_argument} ${synthetic_data_argument}
    """
}
