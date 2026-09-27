process SPLIT_CELL_TYPES {
    tag 'split merged single-cell H5AD'
    publishDir "${params.outdir}/cell_type_splits", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.endsWith('.h5ad') ? filename.tokenize('/').last() : null }

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

process RENDER_CELL_TYPE_REPORT {
    tag { slug }
    publishDir "${params.outdir}/cell_type_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(slug), path(primitive), path(analysis_dir, stageAs: 'analysis_result')
    path pipeline_config
    path report_template

    output:
    path "${slug}"

    script:
    """
    mkdir -p ${slug}
    cp -a analysis_result/. ${slug}/
    pbmc-cell-type-report --input ${primitive} --analysis ${slug}/analysis.h5ad \\
      --output-dir ${slug} --template ${report_template} --config ${pipeline_config}
    """
}

process WRITE_CELL_TYPE_MANIFEST {
    tag 'cell-type analysis manifest'
    publishDir "${params.outdir}/cell_type_analysis", mode: 'copy', overwrite: true

    input:
    path report_dirs
    path merged_input
    path pipeline_config
    path models

    output:
    path 'cell_type_manifest.json'

    script:
    """
    pbmc-cell-type-manifest --input ${merged_input} --config ${pipeline_config} --model ${models} \\
      --report-dir ${report_dirs} --output cell_type_manifest.json
    """
}
