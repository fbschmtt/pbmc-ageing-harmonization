process GENERATE_SYNTHETIC_PSEUDOBULK {
    publishDir "${params.outdir}/analysis/combined/synthetic_de", mode: 'copy', overwrite: true
    input:
    path pipeline_config

    output:
    path 'pseudobulk_merged.h5ad', emit: pseudobulk


    script:
    """
    python -m pbmc_pipeline.synthetic_de --config ${pipeline_config} --output pseudobulk_merged.h5ad
    """
}
