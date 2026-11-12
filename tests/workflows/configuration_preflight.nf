nextflow.enable.dsl = 2

include { VALIDATE_CONFIGURATION } from '../../modules/configuration'
params.preflight_config = null

process ACCEPT_VALIDATED_CONFIG {
    publishDir params.outdir, mode: 'copy'
    input:
    path config
    output:
    path 'accepted.json'
    script:
    """
    cp ${config} accepted.json
    """
}
workflow {
    VALIDATE_CONFIGURATION(channel.value(file(params.preflight_config, checkIfExists: true)),
        ['processing', 'cell_type_analysis', 'differential_expression'], [:], [])
    ACCEPT_VALIDATED_CONFIG(VALIDATE_CONFIGURATION.out.config)
}
