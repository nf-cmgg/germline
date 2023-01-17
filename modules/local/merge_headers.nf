process MERGE_HEADERS {
    tag "$meta.id"
    label 'process_single'

    container 'groovy:3.0.14-jdk8'

    input:
    tuple val(meta), path(vcf), path(ped_vcf)

    output:
    tuple val(meta), path("*.header.txt")   , emit: header
    path  "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    merge_headers.groovy --prefix ${prefix} --vcf ${vcf} --ped_vcf ${ped_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        groovy: \$(echo \$(groovy --version 2>&1) | sed 's/^Groovy Version: //; s/ JVM:.*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.header.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        groovy: \$(echo \$(groovy --version 2>&1) | sed 's/^Groovy Version: //; s/ JVM:.*\$//' )
    END_VERSIONS
    """
}
