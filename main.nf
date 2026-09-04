process CONVERT_RDS {
    tag { study_id }

    input:
    tuple val(study_id), path(source), val(assay), path(converter)

    output:
    tuple val(study_id), path("${study_id}.h5ad"), emit: converted

    script:
    """
    Rscript ${converter} \
      --input ${source} \
      --output ${study_id}.h5ad \
      --assay ${assay}
    """
}

process HARMONIZE {
    tag { study_id }
    publishDir "${params.outdir}/harmonized", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.h5ad') ? name : null }
    publishDir "${params.outdir}/reports", mode: 'copy', overwrite: true,
        saveAs: { name -> name.endsWith('.json') ? name : null }

    input:
    tuple val(study_id), path(expression)
    path tracked_dependencies

    output:
    tuple val(study_id), path("${study_id}.h5ad"), path("${study_id}.json"), emit: harmonized

    script:
    """
    pbmc-harmonize \
      --project-root ${params.project_dir} \
      --config config/pipeline.json \
      --study ${study_id} \
      --input ${expression} \
      --output ${study_id}.h5ad \
      --report-output ${study_id}.json
    """
}

process QC_REPORT {
    tag { study_id }
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true

    input:
    tuple val(study_id), path(expression), path(run_report)
    path qc_template

    output:
    path("${study_id}")

    script:
    """
    mkdir -p ${study_id}
    pbmc-qc \
      --project-root ${params.project_dir} \
      --input ${expression} \
      --run-report ${run_report} \
      --output-dir ${study_id} \
      --study ${study_id} \
      --template ${qc_template}
    """
}

workflow {
    def root = file(params.project_dir)
    def pipeline = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
    def document = new groovy.json.JsonSlurper().parse(
        root.resolve(pipeline.studies_config as String)
    )
    def requested = params.studies == 'all' \
        ? document.studies.keySet() as List \
        : params.studies.toString().split(',').collect { value -> value.trim() }
            .findAll { value -> value }
    def unknown = requested.findAll { study_id -> !document.studies.containsKey(study_id) }
    if (unknown) {
        error "Unknown studies: ${unknown.join(', ')}"
    }

    def direct_inputs = []
    def conversion_inputs = []
    def dependencies = [
        root.resolve('config/pipeline.json'),
        root.resolve(pipeline.studies_config as String),
        root.resolve(pipeline.obs_schema as String),
        root.resolve('reports/qc_report.ipynb')
    ]

    requested.each { study_id ->
        def study = document.studies[study_id]
        def relative_input = params.test \
            ? "${pipeline.test_input_root}/${study.test_input}" \
            : study.input as String
        def conversion_source = params.test
            ? study.conversion?.test_source
            : study.conversion?.source
        if (conversion_source) {
            conversion_inputs << tuple(
                study_id,
                file(root.resolve(conversion_source as String), checkIfExists: true),
                study.conversion.assay ?: 'RNA',
                file(root.resolve('scripts/convert_rds.R'), checkIfExists: true)
            )
        } else {
            direct_inputs << tuple(
                study_id,
                file(root.resolve(relative_input), checkIfExists: true)
            )
        }
        study.joins?.each { join -> dependencies << root.resolve(join.path as String) }
        if (study.annotation.method == 'celltypist') {
            study.annotation.levels.each { level ->
                dependencies << root.resolve(pipeline.models["aifi_${level}"] as String)
            }
        }
    }

    def missing = dependencies.unique().findAll { dependency -> !dependency.exists() }
    if (missing) {
        error "Missing pipeline dependencies:\n  ${missing.join('\n  ')}"
    }

    direct_ch = channel.fromList(direct_inputs)
    conversion_ch = channel.fromList(conversion_inputs)
    converted_ch = CONVERT_RDS(conversion_ch)
    expression_ch = direct_ch.mix(converted_ch)
    dependency_ch = channel.value(dependencies.unique().collect { dependency ->
        file(dependency, checkIfExists: true)
    })

    HARMONIZE(expression_ch, dependency_ch)
    if (!params.skip_qc) {
        qc_template_ch = channel.value(file(
            root.resolve('reports/qc_report.ipynb'), checkIfExists: true
        ))
        QC_REPORT(HARMONIZE.out.harmonized, qc_template_ch)
    }
}
