nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION; clockOptions } from './modules/configuration'

include { AGE_CLOCKS } from './workflows/age_clocks'
include { AGGREGATE_AGE_CLOCKS } from './modules/age_clocks'

workflow {
    def root = file(params.project_dir)
    if (params.aggregate_clocks_only) {
        def input_root = params.age_clock_output ?: "${params.outdir}/analysis/per_type"
        def models = channel.fromPath("${input_root}/*/ageing_clock", type: 'dir', checkIfExists: true)
            .filter { directory -> directory.resolve('age_clock.json').exists() }
            .map { directory -> tuple(directory.parent.name, directory) }
        AGGREGATE_AGE_CLOCKS(
            models.toList().map { records -> tuple(records.collect { record -> record[0] }, records.collect { record -> record[1] }) },
            channel.value([]), [], true,
            channel.value(file(root.resolve('scripts/clock_task.py'), checkIfExists: true)),
            channel.value(file(root.resolve('reports/age_clock_cross_type_report.py'), checkIfExists: true)),
        )
    } else {
        VALIDATE_CONFIGURATION(channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true)),
            ['age_clock'], clockOptions(params), [])
        AGE_CLOCKS(
            channel.value(file(params.pseudobulk_input, checkIfExists: true)),
            VALIDATE_CONFIGURATION.out.config,
            params.clock_cell_type ?: '', root,
        )
    }

}
