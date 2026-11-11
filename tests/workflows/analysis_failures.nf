nextflow.enable.dsl = 2

include { DIFFERENTIAL_EXPRESSION } from '../../workflows/differential_expression'

include { GENE_PROFILES } from '../../workflows/gene_profiles'
include { AGE_CLOCKS } from '../../workflows/age_clocks'
include { CELL_TYPE_ANALYSIS } from '../../workflows/cell_type_analysis'

workflow {
    def root = file(params.project_dir)
    def fixture = channel.value(file("${projectDir}/fixtures/pseudobulk.txt", checkIfExists: true))
    def config = channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true))
    DIFFERENTIAL_EXPRESSION(
        fixture, config,
        channel.value(file(root.resolve('reports/differential_expression_report.py'), checkIfExists: true)),
        channel.value(file(root.resolve('reports/trajectory_analysis_report.py'), checkIfExists: true)),
        channel.value(file("${projectDir}/fixtures", checkIfExists: true)), false, true, true,
    )
    GENE_PROFILES(
        fixture, config,
        channel.of(tuple('test-type', 'Test type', 'TEST_GENE', 'test-gene'), tuple('failed-type', 'Failed type', 'TEST_GENE', 'test-gene')),
        channel.value(file("${projectDir}/fixtures", checkIfExists: true)),
        channel.value(file(root.resolve('scripts/profile_gene.py'), checkIfExists: true)),
        channel.value(file(root.resolve('reports/gene_profile_report.py'), checkIfExists: true)),
    )
    AGE_CLOCKS(fixture, config, '', root)
    CELL_TYPE_ANALYSIS(
        fixture, config,
        channel.value(file(root.resolve('reports/cell_type_report.py'), checkIfExists: true)),
        channel.value([file("${projectDir}/fixtures/model.txt", checkIfExists: true)]), [], channel.empty(),
    )
}
