nextflow.enable.dsl = 2

include { LIST_CLOCK_CELL_TYPES; FIT_AGE_CLOCK; RENDER_AGE_CLOCK; AGGREGATE_AGE_CLOCKS } from '../modules/age_clocks'

workflow AGE_CLOCKS {
    take:
    pseudobulk_input
    pipeline_config
    selection
    root

    main:
    LIST_CLOCK_CELL_TYPES(
        pseudobulk_input, pipeline_config, selection,
        channel.value(file(root.resolve('scripts/list_clock_cell_types.py'), checkIfExists: true)),
    )
    def types = LIST_CLOCK_CELL_TYPES.out.types.splitText().map { line ->
        def fields = line.trim().split('\t', 2)
        tuple(fields[0], fields[1])
    }
    def clock_script = channel.value(file(root.resolve('scripts/clock_task.py'), checkIfExists: true))
    FIT_AGE_CLOCK(pseudobulk_input, pipeline_config, types, pseudobulk_input.map { input -> input.toString() }, clock_script)
    RENDER_AGE_CLOCK(
        FIT_AGE_CLOCK.out.model,
        channel.value(file(root.resolve('reports/age_clock_report.py'), checkIfExists: true)),
        channel.value(file(root.resolve('reports/age_clock_comparison_report.py'), checkIfExists: true)),
    )
    AGGREGATE_AGE_CLOCKS(
        FIT_AGE_CLOCK.out.model.toList().map { records -> tuple(records.collect { record -> record[0] }, records.collect { record -> record[1] }) },
        FIT_AGE_CLOCK.out.status.toList(),
        LIST_CLOCK_CELL_TYPES.out.expected.map { expected -> new groovy.json.JsonSlurper().parse(expected) },
        false, clock_script,
        channel.value(file(root.resolve('reports/age_clock_cross_type_report.py'), checkIfExists: true)),
    )

    emit:
    models = FIT_AGE_CLOCK.out.model
    reports = RENDER_AGE_CLOCK.out.reports
    summary = AGGREGATE_AGE_CLOCKS.out.summary
}
