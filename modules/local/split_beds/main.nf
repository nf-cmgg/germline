process SPLIT_BEDS {
    tag "$meta.id"
    label 'process_single'

    conda "anaconda::gawk=5.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.1.0' :
        'quay.io/biocontainers/gawk:5.1.0' }"

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.bed")  , emit: beds
    path('versions.yml')            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    """
    awk -vFS="\t" '{
        if (\$0 ~ /^[^#].*\$/) {
            if (\$1 ~ /^(chr)?[0-9XY]{1,2}\$/) {
                print \$0 > sprintf("${prefix}_%s_%s_%s.bed", \$1, \$2, \$3)
            } else {
                print \$0 > "${prefix}_others.bed"
            }
        }
    }' ${bed}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    //No stub because this process is so small (would be harder to stub than to just let it run)
}