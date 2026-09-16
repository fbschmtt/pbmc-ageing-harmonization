process CONVERT_RDS {
    tag { study_id }

    input:
    tuple val(study_id), path(source), val(assay), path(converter)

    output:
    tuple val(study_id), path("${study_id}.converted.h5ad"), emit: converted

    script:
    """
    Rscript ${converter} --input ${source} --output ${study_id}.converted.h5ad --assay ${assay}
    """
}
