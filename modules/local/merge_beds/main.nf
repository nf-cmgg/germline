process MERGE_BEDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bed, stageAs: "?/*")
    tuple val(meta2), path(fai)

    output:
    tuple val(meta), path('*.bed'), emit: bed
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for FILE in */*.bed.gz;
    do
        if [[ \$FILE != '*/*.bed.gz' ]]
        then
            gzip -d -f \$FILE
        fi
    done;

    awk '{print \$1"\\t"\$2"\\t"\$3 }' */*.bed | bedtools sort -faidx ${fai} | bedtools merge ${args} > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
