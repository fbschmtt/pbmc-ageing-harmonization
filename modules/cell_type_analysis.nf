process SPLIT_CELL_TYPES {
    tag 'split merged single-cell H5AD'
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename ->
            filename.endsWith('.h5ad')
                ? "${filename.tokenize('/').last().replace('.h5ad', '')}/split.h5ad"
                : null
        }

    input:
    path merged_input
    path pipeline_config

    output:
    path 'cell_types/*.h5ad', emit: cells
    path 'cell_types.json', emit: manifest

    script:
    """
    pbmc-cell-type split --input ${merged_input} --output-dir cell_types --manifest cell_types.json --config ${pipeline_config}
    """

    stub:
    """
    mkdir cell_types
    touch cell_types/test-type.h5ad cell_types/failed-type.h5ad
    printf '%s\n' '{"cell_types":[{"slug":"test-type","cell_type":"Test type"},{"slug":"failed-type","cell_type":"Failed type"}]}' > cell_types.json
    """
}

process ANALYSE_CELL_TYPE {
    tag { slug }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${slug}/analysis/${filename.tokenize('/').last()}" }

    input:
    tuple val(slug), path(primitive)
    path pipeline_config

    output:
    tuple val(slug), path("${slug}/*"), emit: analysis

    script:
    """
    mkdir -p ${slug}
    pbmc-cell-type analyse --input ${primitive} --output-dir ${slug} --config ${pipeline_config}
    """

    stub:
    def fail = params.stub_fail_all_cell_types || (params.stub_fail_cell_type && slug == 'failed-type') ? 'exit 23' : ''
    """
    ${fail}
    mkdir -p ${slug}
    touch ${slug}/analysis.h5ad
    printf '%s\n' '{"status":"complete","cell_type":"Test type","aifi_l1_parent":"Monocyte","n_cells":50,"n_clusters":2}' > ${slug}/report.json
    """
}

workflow ANALYSE_EXISTING_CELL_TYPE_SPLITS {
    take:
    primitives
    pipeline_config

    main:
    ANALYSE_CELL_TYPE(primitives, pipeline_config)

    emit:
    primitives.join(ANALYSE_CELL_TYPE.out.analysis)
}

workflow RENDER_CELL_TYPE_REPORTS {
    take:
    report_inputs
    pipeline_config
    report_template

    main:
    RENDER_CELL_TYPE_REPORT(report_inputs, pipeline_config, report_template)

    emit:
    RENDER_CELL_TYPE_REPORT.out
}

process RENDER_CELL_TYPE_REPORT {
    tag { slug }
    errorStrategy 'ignore'
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${filename}/analysis" }

    input:
    tuple val(slug), path(primitive), path(analysis_files, stageAs: 'analysis_result/*'),
        val(has_merged_de), path(merged_de_input, stageAs: 'merged_de_input')
    path pipeline_config
    path report_template

    output:
    path "${slug}"

    script:
    """
    mkdir -p ${slug}
    cp -aL analysis_result/. ${slug}/
    if [ "${has_merged_de}" = "true" ]; then
      mkdir -p differential_expression
      cp merged_de_input differential_expression/merged.csv
    fi
    pbmc-cell-type-report --input ${primitive} --analysis ${slug}/analysis.h5ad \\
      --output-dir ${slug} --template ${report_template} --config ${pipeline_config}
    """

    stub:
    def fail = params.stub_fail_cell_report && slug == 'failed-type' ? 'exit 41' : ''
    """
    ${fail}
    mkdir -p ${slug}
    cp -aL analysis_result/. ${slug}/
    if [ "${has_merged_de}" = "true" ]; then
      cp merged_de_input ${slug}/received_merged_de.csv
    fi
    printf '%s\n' '<html>Cell report stub</html>' > ${slug}/report.html
    printf '%s\n' '{}' > ${slug}/executed.ipynb
    """
}

process WRITE_CELL_TYPE_MANIFEST {
    errorStrategy 'ignore'
    tag 'cell-type analysis manifest'
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true

    input:
    path report_dirs
    path split_manifest
    val analysed_slugs
    val selected_slugs
    path merged_input
    path pipeline_config
    path models

    output:
    path 'cell_type_manifest.json'

    script:
    def analysed_args = analysed_slugs
        ? analysed_slugs.collect { slug -> "--analysed-slug ${slug}" }.join(' ')
        : ''
    def selected_args = selected_slugs
        ? selected_slugs.collect { slug -> "--selected-slug ${slug}" }.join(' ')
        : ''
    """
    pbmc-cell-type-manifest --input ${merged_input} --config ${pipeline_config} --model ${models} \\
      --split-manifest ${split_manifest} ${analysed_args} ${selected_args} \\
      --report-dir ${report_dirs} --output cell_type_manifest.json
    """

    stub:
    def analysed_args = analysed_slugs.collect { slug -> "--analysed-slug ${slug}" }.join(' ')
    def selected_args = selected_slugs.collect { slug -> "--selected-slug ${slug}" }.join(' ')
    """
    python3 ${params.project_dir}/src/pbmc_pipeline/cell_type_manifest.py \
      --input ${merged_input} --config ${pipeline_config} --model ${models} \
      --split-manifest ${split_manifest} ${analysed_args} ${selected_args} \
      --report-dir ${report_dirs} --output cell_type_manifest.json
    """
}
