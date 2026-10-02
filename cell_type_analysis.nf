nextflow.enable.dsl = 2

include {
    SPLIT_CELL_TYPES
    ANALYSE_CELL_TYPE
    RENDER_CELL_TYPE_REPORT
    WRITE_CELL_TYPE_MANIFEST
    SPLIT_AND_ANALYSE_CELL_TYPES
    ANALYSE_EXISTING_CELL_TYPE_SPLITS
    RENDER_CELL_TYPE_REPORTS
} from './modules/cell_type_analysis'

workflow {
    def root = file(params.project_dir)
    def pipeline_config_file = file(root.resolve('config/pipeline.json'), checkIfExists: true)
    def pipeline_config = channel.value(pipeline_config_file)
    def report_template = channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true))
    def de_manifest_path = params.differential_expression_manifest
        ?: "${params.outdir}/differential_expression/differential_expression.json"
    def de_manifest_file = file(de_manifest_path.toString())
    def de_by_slug = channel.empty()
    if (de_manifest_file.isFile()) {
        def de_manifest = new groovy.json.JsonSlurper().parse(de_manifest_file)
        def result_index = de_manifest.results_by_cell_type
        if (result_index == null && de_manifest.cell_types instanceof List) {
            // Read manifests produced before the explicit artifact index existed.
            result_index = de_manifest.cell_types.collectEntries { record ->
                def completed = (record.models ?: []).findAll { model -> model.status == 'complete' }
                [(record.slug as String): [
                    per_study: completed.findAll { model -> model.model == 'per_study' }
                        .collect { model -> [study: model.study, path: model.results] },
                    merged: completed.find { model -> model.model == 'merged' }?.results,
                ]]
            }
        }
        if (!(result_index instanceof Map)) {
            error "${de_manifest_file}: missing results_by_cell_type DE artifact index"
        }
        def indexed_dirs = result_index.collect { slug, record ->
            def result_paths = (record.per_study ?: []).collectMany { model ->
                if (model.paths instanceof Map) {
                    return model.paths.values().collect { result_path -> result_path as String }
                }
                return model.path ? [model.path as String] : []
            } + (record.combined instanceof Map
                ? record.combined.values().collect { result_path -> result_path as String }
                : []) + (record.merged ? [record.merged as String] : [])
            result_paths = result_paths.unique()
            result_paths.each { relative_path ->
                if (!de_manifest_file.parent.resolve(relative_path as String).isFile()) {
                    error "${de_manifest_file}: indexed DE result is missing: ${relative_path}"
                }
            }
            tuple(slug as String, file(de_manifest_file.parent.resolve(slug as String), checkIfExists: true))
        }
        de_by_slug = channel.fromList(indexed_dirs)
    }

    def report_inputs = null
    def merged_input = null
    def models = null
    if (params.cell_type_analysis_input) {
        if (!params.cell_type_input) {
            error '--cell_type_analysis_input requires --cell_type_input'
        }
        if (params.split_only) {
            error '--cell_type_analysis_input and --split_only cannot be used together'
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.baseName, primitive) }
        def analyses = channel.fromPath(
            params.cell_type_analysis_input, type: 'dir', checkIfExists: true,
        ).map { analysis_dir -> tuple(analysis_dir.baseName, analysis_dir) }
        report_inputs = primitives.join(analyses)
    } else if (params.cell_type_input) {
        if (params.split_only) {
            error '--cell_type_input and --split_only cannot be used together'
        }
        def primitives = channel.fromPath(params.cell_type_input, checkIfExists: true)
            .map { primitive -> tuple(primitive.baseName, primitive) }
        ANALYSE_EXISTING_CELL_TYPE_SPLITS(primitives, pipeline_config)
        report_inputs = ANALYSE_EXISTING_CELL_TYPE_SPLITS.out
    } else {
        if (!params.merged_input) {
            error 'Provide --merged_input, for example output/merged/single_cell_merged.h5ad'
        }
        def pipeline = new groovy.json.JsonSlurper().parse(pipeline_config_file)
        merged_input = channel.value(file(params.merged_input, checkIfExists: true))
        if (params.split_only) {
            SPLIT_CELL_TYPES(merged_input, pipeline_config)
        } else {
            models = channel.value([
                file(root.resolve(pipeline.models.aifi_l1 as String), checkIfExists: true),
                file(root.resolve(pipeline.models.aifi_l2 as String), checkIfExists: true),
            ])
            SPLIT_AND_ANALYSE_CELL_TYPES(merged_input, pipeline_config)
            report_inputs = SPLIT_AND_ANALYSE_CELL_TYPES.out
        }
    }

    if (report_inputs != null) {
        RENDER_CELL_TYPE_REPORTS(
            report_inputs,
            de_by_slug,
            pipeline_config_file,
            pipeline_config,
            report_template,
        )
        if (models != null) {
            WRITE_CELL_TYPE_MANIFEST(
                RENDER_CELL_TYPE_REPORTS.out.collect(),
                merged_input,
                pipeline_config,
                models,
            )
        }
    }
}
