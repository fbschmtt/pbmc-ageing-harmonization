process SPLIT_CELL_TYPES {
    tag 'split merged single-cell H5AD'

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

process RUN_CELL_TYPE_REPORT {
    tag { primitive.baseName }
    publishDir "${params.outdir}/cell_type_analysis", mode: 'copy', overwrite: true

    input:
    path primitive
    path pipeline_config
    path report_template

    output:
    path "${primitive.baseName}"

    script:
    """
    mkdir -p ${primitive.baseName}
    pbmc-cell-type-report --input ${primitive} --output-dir ${primitive.baseName} --template ${report_template} --config ${pipeline_config}
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
