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

    stub:
    """
    printf 'test-type\tTest type\nfailed-type\tFailed type\n' > cell_types.tsv
    """
}

process RUN_CELL_TYPE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION {
    errorStrategy 'ignore'
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
    tuple val(slug), path("cell_type_result/${slug}"), emit: results

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir cell_type_result --config ${pipeline_config} --cell-type "${cell_type}" \
      --cpus ${task.cpus} ${test_mode_argument} ${synthetic_data_argument}
    """

    stub:
    def fail = params.stub_fail_de && slug == 'failed-type' ? 'exit 47' : ''
    """
    ${fail}
    mkdir -p cell_type_result/${slug}
    printf '%s\n' '{"cell_type":"${cell_type}","analysis_mode":"test_only"}' > cell_type_result/${slug}/run_metadata.json
    printf 'gene,log2FoldChange\nTEST_GENE,0.25\n' > cell_type_result/${slug}/merged.csv
    """
}

process WRITE_PSEUDOBULK_DIFFERENTIAL_EXPRESSION_MANIFEST {
    errorStrategy 'ignore'
    tag 'collect cell-type PyDESeq2 models'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'analysis/per_type' ? null : filename }

    input:
    path pseudobulk_input
    path pipeline_config
    path cell_type_results
    val de_test_mode
    val synthetic_test_data
    val allow_partial_results

    output:
    path 'analysis/combined/differential_expression', emit: results
    path 'analysis/per_type', emit: per_type

    script:
    def test_mode_argument = de_test_mode ? '--test-mode' : ''
    def synthetic_data_argument = synthetic_test_data ? '--synthetic-test-data' : ''
    def partial_results_argument = allow_partial_results ? '--allow-partial-results' : ''
    """
    pbmc-differential-expression --input ${pseudobulk_input} \
      --output-dir analysis/combined/differential_expression \
      --per-type-output-dir analysis/per_type --config ${pipeline_config} \
      --combine-cell-type-results ${cell_type_results} \
      ${test_mode_argument} ${synthetic_data_argument} ${partial_results_argument}
    """

    stub:
    """
    mkdir -p analysis/combined/differential_expression analysis/per_type
    for result in ${cell_type_results}; do
      slug=\$(basename \"\$result\")
      mkdir -p analysis/per_type/\$slug/differential_expression
      cp -a \"\$result\"/. analysis/per_type/\$slug/differential_expression/
    done
    printf '%s\n' '{"test_only":true}' > analysis/combined/differential_expression/differential_expression.json
    """
}

process RENDER_CROSS_CELL_TYPE_ANALYSIS_REPORT {
    errorStrategy 'ignore'
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

    stub:
    """
    mkdir trajectory
    printf '%s\n' '<html>Cross-type stub</html>' > trajectory/report.html
    """
}

process ANALYSE_CELL_TYPE_DE_DIAGNOSTICS {
    errorStrategy 'ignore'
    tag { slug }
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            filename.startsWith('diagnostics/')
                ? "${slug}/differential_expression/${filename.substring('diagnostics/'.length())}"
                : null
        }

    input:
    tuple val(slug), path(result_dir, stageAs: 'result')
    path pipeline_config
    path go_enrichment_data

    output:
    tuple val(slug), path('diagnostics/*'), emit: diagnostics

    script:
    """
    PBMC_GO_ENRICHMENT_DATA_DIR="\$PWD/${go_enrichment_data}" pbmc-de-diagnostics \\
      --result-dir result --output-dir diagnostics --config ${pipeline_config}
    """

    stub:
    """
    test -f result/run_metadata.json
    mkdir -p diagnostics/residual_wgcna
    printf '%s\n' '{"schema_version":1,"test_only":true}' > diagnostics/diagnostics.json
    printf 'study,estimate\nTEST_STUDY,0.1\n' > diagnostics/fraction_model_coefficients.tsv
    printf '%s\n' '{"status":"test_only"}' > diagnostics/residual_wgcna/metadata.json
    """
}

process RENDER_PER_TYPE_DIFFERENTIAL_EXPRESSION_REPORTS {
    errorStrategy 'ignore'
    tag { slug }
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            filename.startsWith("${slug}/differential_expression/") ? filename : null
        }

    input:
    tuple val(slug), path(result_dir, stageAs: 'result'), path(diagnostic_files, stageAs: 'diagnostics/*')
    path pipeline_config
    path report_template

    output:
    path "${slug}/differential_expression/{executed.ipynb,report.html,report_render.log}", emit: reports

    script:
    """
    mkdir -p ${slug}/differential_expression
    pbmc-differential-expression-report --result-dir result --analysis-dir diagnostics \\
      --output-dir ${slug}/differential_expression \\
      --template ${report_template} --config ${pipeline_config} --project-root .
    """

    stub:
    def fail = params.stub_fail_de_report && slug == 'failed-type' ? 'exit 43' : ''
    """
    ${fail}
    test -f result/run_metadata.json
    test -f diagnostics/diagnostics.json
    test -f diagnostics/residual_wgcna/metadata.json
    mkdir -p ${slug}/differential_expression
    printf '%s\n' '<html>Publication test only</html>' > ${slug}/differential_expression/report.html
    printf '%s\n' '{"cells":[]}' > ${slug}/differential_expression/executed.ipynb
    printf '%s\n' 'publication test only' > ${slug}/differential_expression/report_render.log
    """
}
