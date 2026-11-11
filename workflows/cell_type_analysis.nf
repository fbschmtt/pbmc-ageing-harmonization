nextflow.enable.dsl = 2

include { SPLIT_CELL_TYPES; ANALYSE_CELL_TYPE; RENDER_CELL_TYPE_REPORTS; WRITE_CELL_TYPE_MANIFEST } from '../modules/cell_type_analysis'
include { selectCellTypes } from '../modules/cell_type_selection'

workflow CELL_TYPE_ANALYSIS {
    take:
    merged_input
    pipeline_config
    report_template
    models
    selected_slugs
    de_results

    main:
    SPLIT_CELL_TYPES(merged_input, pipeline_config)
    def primitives = SPLIT_CELL_TYPES.out.cells.flatten()
        .map { primitive -> tuple(primitive.baseName, primitive) }
        .toList()
        .flatMap { entries ->
            def selection = selectCellTypes(entries, selected_slugs)
            if (selection.error) error(selection.error)
            selection.entries
        }
    ANALYSE_CELL_TYPE(primitives, pipeline_config)
    // The barrier covers only current-run producers. No lookup in the published tree.
    def current_de = de_results.toList().map { records ->
        records.collectEntries { slug, directory -> [(slug.toString()): directory.resolve('merged.csv')] }
    }
    def reports_with_config = primitives.join(ANALYSE_CELL_TYPE.out.analysis)
        .combine(current_de).combine(pipeline_config)
        .map { slug, primitive, analysis_dir, de, config_file ->
            def merged_de = de[slug.toString()]
            def available = merged_de != null && merged_de.exists()
            tuple(slug, primitive, analysis_dir, available, available ? merged_de : config_file)
        }
    RENDER_CELL_TYPE_REPORTS(reports_with_config, pipeline_config, report_template)
    WRITE_CELL_TYPE_MANIFEST(
        RENDER_CELL_TYPE_REPORTS.out.toList(), SPLIT_CELL_TYPES.out.manifest,
        ANALYSE_CELL_TYPE.out.analysis.map { slug, _directory -> slug.toString() }.toList(),
        selected_slugs, merged_input, pipeline_config, models,
    )

    emit:
    reports = RENDER_CELL_TYPE_REPORTS.out
    manifest = WRITE_CELL_TYPE_MANIFEST.out
}
