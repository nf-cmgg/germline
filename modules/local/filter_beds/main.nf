process FILTER_BEDS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'quay.io/biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path('*.bed.gz'), emit: bed
    path  "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Remove regions with no coverage from the callable regions BED file
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat ${bed} | grep -v NO_COVERAGE | bgzip --threads ${task.cpus} --stdout > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
