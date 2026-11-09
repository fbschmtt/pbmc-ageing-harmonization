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
}

process ANALYSE_CELL_TYPE {
    tag { slug }
    errorStrategy { params.merged_input && !params.split_only ? 'ignore' : 'terminate' }

    input:
    tuple val(slug), path(primitive)
    path pipeline_config

    output:
    tuple val(slug), path("${slug}"), emit: analysis

    script:
    """
    mkdir -p ${slug}
    pbmc-cell-type analyse --input ${primitive} --output-dir ${slug} --config ${pipeline_config}
    """
}

workflow SPLIT_AND_ANALYSE_CELL_TYPES {
    take:
    merged_input
    pipeline_config

    main:
    SPLIT_CELL_TYPES(merged_input, pipeline_config)
    def primitives = SPLIT_CELL_TYPES.out.cells.flatten()
        .map { primitive -> tuple(primitive.baseName, primitive) }
    ANALYSE_CELL_TYPE(primitives, pipeline_config)

    emit:
    primitives.join(ANALYSE_CELL_TYPE.out.analysis)
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
    errorStrategy { params.merged_input && !params.split_only ? 'ignore' : 'terminate' }
    publishDir "${params.outdir}/analysis/per_type", mode: 'copy', overwrite: true,
        saveAs: { filename -> "${filename}/analysis" }

    input:
    tuple val(slug), path(primitive), path(analysis_dir, stageAs: 'analysis_result'),
        val(has_merged_de), path(merged_de_input, stageAs: 'merged_de_input')
    path pipeline_config
    path report_template

    output:
    path "${slug}"

    script:
    """
    mkdir -p ${slug}
    cp -a analysis_result/. ${slug}/
    if [ "${has_merged_de}" = "true" ]; then
      mkdir -p differential_expression
      cp merged_de_input differential_expression/merged.csv
    fi
    pbmc-cell-type-report --input ${primitive} --analysis ${slug}/analysis.h5ad \\
      --output-dir ${slug} --template ${report_template} --config ${pipeline_config}
    """
}

process WRITE_CELL_TYPE_MANIFEST {
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
}
