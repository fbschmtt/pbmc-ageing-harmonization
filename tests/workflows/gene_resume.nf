nextflow.enable.dsl = 2

include { GENE_PROFILES } from '../../workflows/gene_profiles'

workflow {
    def root = file(params.project_dir)
    GENE_PROFILES(
        channel.value(file("${projectDir}/fixtures/pseudobulk.txt", checkIfExists: true)),
        channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true)),
        channel.of(tuple('test-type', 'Test type', 'GENE_1', 'gene-1'), tuple('test-type', 'Test type', 'GENE_2', 'gene-2')),
        channel.value(file("${projectDir}/fixtures", checkIfExists: true)),
        channel.value(file(root.resolve('scripts/profile_gene.py'), checkIfExists: true)),
        channel.value(file(params.resume_template, checkIfExists: true)),
    )
}
