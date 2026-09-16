nextflow.enable.dsl = 2

include { CONVERT_RDS } from './modules/conversion'
include { HARMONIZE } from './modules/harmonization'
include { PSEUDOBULK; MERGE_PSEUDOBULKS; MERGE_SINGLE_CELLS } from './modules/merge'
include { RENDER_STUDY_QC; RENDER_MERGE_QC } from './modules/qc'
include { WRITE_RUN_MANIFEST } from './modules/manifest'

workflow {
    def root = file(params.project_dir)
    def pipeline = new groovy.json.JsonSlurper().parse(root.resolve('config/pipeline.json'))
    def document = new groovy.json.JsonSlurper().parse(root.resolve(pipeline.studies_config as String))
    def requested = params.studies == 'all'
        ? document.studies.keySet() as List
        : params.studies.toString().split(',').collect { value -> value.trim() }.findAll { value -> value }
    def unknown = requested.findAll { study_id -> !document.studies.containsKey(study_id) }
    if (unknown) {
        error "Unknown studies: ${unknown.join(', ')}"
    }

    def direct_inputs = []
    def conversion_inputs = []
    def dependencies = [
        root.resolve('config/pipeline.json'), root.resolve(pipeline.studies_config as String),
        root.resolve(pipeline.obs_schema as String), root.resolve('reports/qc_report.ipynb'),
        root.resolve('reports/merge_qc_report.ipynb'), root.resolve(pipeline.models.aifi_l2 as String)
    ]
    requested.each { study_id ->
        def study = document.studies[study_id]
        def input = params.test ? "${pipeline.test_input_root}/${study.test_input}" : study.input as String
        def conversion_source = params.test ? study.conversion?.test_source : study.conversion?.source
        if (conversion_source) {
            conversion_inputs << tuple(study_id, file(root.resolve(conversion_source as String), checkIfExists: true),
                study.conversion.assay ?: 'RNA', file(root.resolve('scripts/convert_rds.R'), checkIfExists: true))
        } else {
            direct_inputs << tuple(study_id, file(root.resolve(input), checkIfExists: true))
        }
        study.joins?.each { join -> dependencies << root.resolve(join.path as String) }
        if (study.annotation.method == 'celltypist') {
            study.annotation.levels.each { level -> dependencies << root.resolve(pipeline.models["aifi_${level}"] as String) }
        }
    }
    def missing = dependencies.unique().findAll { dependency -> !dependency.exists() }
    if (missing) {
        error "Missing pipeline dependencies:\n  ${missing.join('\n  ')}"
    }

    def dependency_ch = channel.value(dependencies.unique().collect { dependency -> file(dependency, checkIfExists: true) })
    def merge_config_ch = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    def aifi_l2_model_ch = channel.value(file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true))
    def expression_ch = channel.fromList(direct_inputs).mix(CONVERT_RDS(channel.fromList(conversion_inputs)))

    HARMONIZE(expression_ch, dependency_ch)
    PSEUDOBULK(HARMONIZE.out.harmonized, merge_config_ch)
    MERGE_PSEUDOBULKS(PSEUDOBULK.out.pseudobulk.map { _study_id, expression, _report -> expression }.collect(), merge_config_ch)

    def merged_outputs = MERGE_PSEUDOBULKS.out.merged
    if (params.merge_single_cell || pipeline.merge.single_cell.enabled) {
        MERGE_SINGLE_CELLS(HARMONIZE.out.harmonized.map { _study_id, expression, _report -> expression }.collect(), merge_config_ch, aifi_l2_model_ch)
        merged_outputs = merged_outputs.mix(MERGE_SINGLE_CELLS.out.merged)
    }
    WRITE_RUN_MANIFEST(
        merged_outputs.map { _name, expression, _report -> expression }.collect(),
        merged_outputs.map { _name, _expression, report -> report }.collect(),
        merge_config_ch,
        requested.join(','),
    )
    if (!params.skip_qc) {
        def study_qc_template = channel.value(file(root.resolve('reports/qc_report.ipynb'), checkIfExists: true))
        RENDER_STUDY_QC(HARMONIZE.out.harmonized, study_qc_template)
        def merge_qc_template = channel.value(file(root.resolve('reports/merge_qc_report.ipynb'), checkIfExists: true))
        RENDER_MERGE_QC(
            merged_outputs.map { _name, expression, _report -> expression }.collect(),
            merged_outputs.map { _name, _expression, report -> report }.collect(),
            merge_qc_template,
        )
    }
}
