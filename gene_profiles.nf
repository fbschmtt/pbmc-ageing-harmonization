nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION; profileCountError } from './modules/configuration'

include { TOP_GENE_PROFILES; GENE_PROFILES } from './workflows/gene_profiles'

workflow {
    if (!params.profile_gene) {
        def problem = profileCountError(params.top_age_gene_count)
        if (problem) error(problem)
    }
    def root = file(params.project_dir)
    VALIDATE_CONFIGURATION(channel.value(file(root.resolve('config/pipeline.json'), checkIfExists: true)), ['differential_expression'], [:], [])
    def config = VALIDATE_CONFIGURATION.out.config
    def input = channel.value(file(params.pseudobulk_input, checkIfExists: true))
    def de_root = file(params.per_type_de_dir ?: "${params.outdir}/analysis/per_type", checkIfExists: false)
    def context = channel.value(de_root.exists() ? de_root : file(root.resolve('config'), checkIfExists: true))
    if (params.profile_gene) {
        if (!params.profile_cell_type) error '--profile_gene requires --profile_cell_type'
        def label = params.profile_cell_type.toString()
        def gene = params.profile_gene.toString()
        def slug = label.toLowerCase().replaceAll(/[^a-z0-9]+/, '-').replaceAll(/^-|-$/, '') ?: 'unnamed'
        def gene_slug = gene.toLowerCase().replaceAll(/[^a-z0-9]+/, '-').replaceAll(/^-|-$/, '') ?: 'unnamed'
        GENE_PROFILES(
            input, config, channel.of(tuple(slug, label, gene, gene_slug)), context,
            channel.value(file(root.resolve('scripts/profile_gene.py'), checkIfExists: true)),
            channel.value(file(root.resolve('reports/gene_profile_report.py'), checkIfExists: true)),
        )
    } else {
        def models = channel.fromPath("${de_root}/*/differential_expression", type: 'dir', checkIfExists: true)
            .map { directory -> tuple(directory.parent.name, directory) }
            .combine(config).map { slug, directory, _config -> tuple(slug, directory) }
        TOP_GENE_PROFILES(input, config, models, params.top_age_gene_count as int, root)
    }
}
