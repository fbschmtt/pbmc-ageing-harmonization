nextflow.enable.dsl = 2

include { COLLECT_PROFILE_DE_CONTEXT; SELECT_PROFILE_GENES; PREPARE_GENE_PROFILE_MODEL; FIT_GENE_PROFILE; RENDER_GENE_PROFILE } from '../modules/gene_profiles'

workflow GENE_PROFILES {
    take:
    pseudobulk_input
    pipeline_config
    requests
    per_type_de
    profile_script
    report_template

    main:
    PREPARE_GENE_PROFILE_MODEL(
        pseudobulk_input, pipeline_config,
        requests.map { slug, cell_type, _gene, _gene_slug -> tuple(slug, cell_type) }.unique(),
        profile_script,
    )
    FIT_GENE_PROFILE(
        pseudobulk_input, pipeline_config,
        requests.combine(PREPARE_GENE_PROFILE_MODEL.out.model, by: 0),
        per_type_de, profile_script,
    )
    RENDER_GENE_PROFILE(FIT_GENE_PROFILE.out.profile, profile_script, report_template)

    emit:
    data = FIT_GENE_PROFILE.out.profile
    reports = RENDER_GENE_PROFILE.out.reports
}

workflow TOP_GENE_PROFILES {
    take:
    pseudobulk_input
    pipeline_config
    de_results
    top_n
    root

    main:
    COLLECT_PROFILE_DE_CONTEXT(de_results.map { _slug, directory -> directory }.toList())
    SELECT_PROFILE_GENES(de_results, top_n, channel.value([
        file(root.resolve('scripts/select_gene_profiles.py'), checkIfExists: true),
        file(root.resolve('scripts/profile_top_age_genes.py'), checkIfExists: true),
    ]))
    def selected_requests = SELECT_PROFILE_GENES.out.selection.flatMap { slug, selection ->
        selection.readLines().findAll { line -> line.trim() }.collect { line ->
            def fields = line.split('\t', -1)
            if (fields.size() != 2) error "Invalid gene selection row: ${line}"
            def gene_slug = fields[1].toLowerCase().replaceAll(/[^a-z0-9]+/, '-').replaceAll(/^-|-$/, '') ?: 'unnamed'
            tuple(slug, fields[0], fields[1], gene_slug)
        }
    }
    def requests = params.verify_gene_profile
        ? selected_requests.mix(channel.of(tuple('cd14-monocyte', 'CD14 monocyte', 'SYNTH_AGE_MARKER_1', 'synth-age-marker-1'))).unique()
        : selected_requests
    GENE_PROFILES(
        pseudobulk_input, pipeline_config, requests, COLLECT_PROFILE_DE_CONTEXT.out.context,
        channel.value(file(root.resolve('scripts/profile_gene.py'), checkIfExists: true)),
        channel.value(file(root.resolve('reports/gene_profile_report.py'), checkIfExists: true)),
    )

    emit:
    data = GENE_PROFILES.out.data
    reports = GENE_PROFILES.out.reports
}
