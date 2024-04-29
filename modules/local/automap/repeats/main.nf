process AUTOMAP_REPEATS {
    tag "$meta.id"
    label 'process_single'

    container "cmgg/automap:1.0.0"

    input:
    tuple val(meta), val(genome)

    output:
    tuple val(meta), path("*.bed")  , emit: repeats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = "1.0.0"

    // Files are present in the container
    if (genome == "hg38" || genome == "GRCh38") {
        """
        zcat /usr/local/lib/automap/Resources/repeats_hg38.part*.bed.gz > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            automap: $VERSION
        END_VERSIONS
        """
    } else if (genome == "hg19" || genome == "GRCh37") {
        """
        zcat /usr/local/lib/automap/Resources/repeats.part*.bed.gz > ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            automap: $VERSION
        END_VERSIONS
        """
    } else {
        error("No repeat regions can be found for genome '${genome}'. Available genomes are hg38/GRCh38 or hg19/GRCh37")
    }


    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = "1.0.0"
    if (["hg38", "GRCh38", "hg19", "GRCh37"].contains(genome)) {
        """
        touch ${prefix}.bed

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            automap: $VERSION
        END_VERSIONS
        """
    } else {
        error("No repeat regions can be found for genome '${genome}'. Available genomes are hg38/GRCh38 or hg19/GRCh37")
    }
}
