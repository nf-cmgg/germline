process MERGE_VCF_HEADERS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.9" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9' :
        'quay.io/biocontainers/python:3.9' }"

    input:
    tuple val(meta), path(vcf1), path(vcf2)

    output:
    tuple val(meta), path("*.vcf")  , emit: vcf
    path "versions.yml"             , emit: versions

    script: // This script is bundled with the pipeline, in TVA/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    merge_vcf_headers.py \\
        $vcf1 \\
        $vcf2 \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
