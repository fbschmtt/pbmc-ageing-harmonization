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
    de_by_slug
    pipeline_config_file
    pipeline_config
    report_template

    main:
    def report_routes = report_inputs.join(de_by_slug, remainder: true).branch { row ->
        with_de: row.size() > 3 && row[1] != null && row[3] != null
        without_de: row.size() > 1 && row[1] != null && (row.size() < 4 || row[3] == null)
    }
    def with_de = report_routes.with_de.map { row ->
        tuple(row[0], row[1], row[2], row[3], true)
    }
    def without_de = report_routes.without_de.map { row ->
        // Reuse a known staged file in the uniform path slot; the flag prevents
        // the report command from treating it as a DE result directory.
        tuple(row[0], row[1], row[2], pipeline_config_file, false)
    }
    RENDER_CELL_TYPE_REPORT(with_de.mix(without_de), pipeline_config, report_template)

    emit:
    RENDER_CELL_TYPE_REPORT.out
}

process RENDER_CELL_TYPE_REPORT {
    tag { slug }
    publishDir "${params.outdir}/cell_type_analysis", mode: 'copy', overwrite: true

    input:
    tuple val(slug), path(primitive), path(analysis_dir, stageAs: 'analysis_result'),
        path(de_payload, stageAs: 'de_input'), val(has_de)
    path pipeline_config
    path report_template

    output:
    path "${slug}"

    script:
    def de_argument = has_de ? "--differential-expression-dir ${de_payload}" : ''
    """
    mkdir -p ${slug}
    cp -a analysis_result/. ${slug}/
    pbmc-cell-type-report --input ${primitive} --analysis ${slug}/analysis.h5ad \\
      --output-dir ${slug} --template ${report_template} --config ${pipeline_config} \\
      ${de_argument}
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
