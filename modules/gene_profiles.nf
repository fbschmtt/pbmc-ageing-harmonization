process SELECT_PROFILE_GENES {
    tag { slug }
    errorStrategy 'ignore'

    input:
    tuple val(slug), path(result_dir, stageAs: 'result')
    val top_n
    path selection_scripts

    output:
    tuple val(slug), path('genes.tsv'), emit: selection

    script:
    """
    python select_gene_profiles.py --result-dir result --top-n ${top_n} --output genes.tsv
    """

    stub:
    """
    printf 'Test type\tTEST_GENE\n' > genes.tsv
    """
}

process PREPARE_GENE_PROFILE_MODEL {
    tag { slug }
    errorStrategy 'ignore'
    maxForks 2

    input:
    path pseudobulk_input
    path pipeline_config
    tuple val(slug), val(cell_type)
    path profile_script

    output:
    tuple val(slug), path('model'), emit: model

    script:
    """
    python ${profile_script} --input ${pseudobulk_input} --config ${pipeline_config} \\
      --cell-type '${cell_type.replace("'", "'\"'\"'")}' --prepare-model --output-dir model --cpus ${task.cpus}
    """

    stub:
    """
    mkdir model
    printf '%s\n' '{"test_only":true}' > model/manifest.json
    """
}

process FIT_GENE_PROFILE {
    tag { "${slug}/${gene_slug}" }
    errorStrategy 'ignore'
    maxForks 2
    publishDir { params.gene_profile_output ?: "${params.outdir}/analysis/per_gene" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.startsWith('profiles/') ? filename.substring('profiles/'.length()) : null }

    input:
    path pseudobulk_input
    path pipeline_config
    tuple val(slug), val(cell_type), val(gene), val(gene_slug), path(model_dir, stageAs: 'model')
    path per_type_de, stageAs: 'current_de'
    path profile_script

    output:
    tuple val(slug), val(cell_type), val(gene_slug), path("profiles/${slug}/${gene_slug}/*"), emit: profile

    script:
    """
    python ${profile_script} --input ${pseudobulk_input} --config ${pipeline_config} \\
      --cell-type '${cell_type.replace("'", "'\"'\"'")}' --gene '${gene.replace("'", "'\"'\"'")}' \\
      --prepared-model model --no-render --output-dir profiles --differential-expression-dir current_de --cpus ${task.cpus}
    """

    stub:
    def fail = params.stub_fail_gene && slug == 'failed-type' ? 'exit 17' : ''
    """
    ${fail}
    mkdir -p profiles/${slug}/${gene_slug}
    printf '%s\n' '{"gene":"${gene}","cell_type":"${cell_type}"}' > profiles/${slug}/${gene_slug}/fit_summary.json
    printf 'gene,estimate\n${gene},0.1\n' > profiles/${slug}/${gene_slug}/gene_age_result.csv
    """
}

process RENDER_GENE_PROFILE {
    tag { "${slug}/${gene_slug}" }
    errorStrategy 'ignore'
    publishDir { params.gene_profile_output ?: "${params.outdir}/analysis/per_gene" }, mode: 'copy', overwrite: true,
        saveAs: { filename -> filename.startsWith('profiles/') ? filename.substring('profiles/'.length()) : null }

    input:
    tuple val(slug), val(cell_type), val(gene_slug), path(profile_files, stageAs: 'prepared/*')
    path profile_script
    path report_template

    output:
    path "profiles/${slug}/${gene_slug}/{report.html,executed.ipynb,report_render.log}", emit: reports

    script:
    """
    mkdir -p profiles/${slug}/${gene_slug}
    cp -aL prepared/. profiles/${slug}/${gene_slug}/
    python ${profile_script} --cell-type '${cell_type.replace("'", "'\"'\"'")}' --render-only \\
      --output-dir profiles/${slug}/${gene_slug} --template ${report_template}
    """

    stub:
    def fail = params.stub_fail_gene_report && slug == 'failed-type' ? 'exit 31' : ''
    """
    ${fail}
    mkdir -p profiles/${slug}/${gene_slug}
    printf '%s\n' '<html>Gene profile stub</html>' > profiles/${slug}/${gene_slug}/report.html
    printf '%s\n' '{}' > profiles/${slug}/${gene_slug}/executed.ipynb
    printf '%s\n' 'stub' > profiles/${slug}/${gene_slug}/report_render.log
    """
}

process COLLECT_PROFILE_DE_CONTEXT {
    errorStrategy 'ignore'
    input:
    path result_dirs, stageAs: 'results/*'

    output:
    path 'current_de', emit: context

    script:
    """
    mkdir current_de
    for result in ${result_dirs}; do
      slug=\$(basename \"\$result\")
      mkdir -p current_de/\$slug/differential_expression
      cp -a \"\$result\"/. current_de/\$slug/differential_expression/
    done
    """

    stub:
    """
    mkdir current_de
    """
}
