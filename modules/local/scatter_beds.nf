process SCATTER_BEDS {
    tag "$meta.id"
    label 'process_single'

    container 'groovy:3.0.14-jdk8'

    input:
    tuple val(meta), path(bed)
    val(scatter_size)

    output:
    tuple val(meta), path("*.bed"), emit: scatter

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    scatter_beds.groovy --prefix ${prefix} --bed ${bed} --size ${scatter_size}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        groovy: \$(echo \$(groovy --version 2>&1) | sed 's/^Groovy Version: //; s/ JVM:.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}_1.bed
    touch ${prefix}_2.bed
    touch ${prefix}_3.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        groovy: \$(echo \$(groovy --version 2>&1) | sed 's/^Groovy Version: //; s/ JVM:.*\$//' )
    END_VERSIONS
    """
}
