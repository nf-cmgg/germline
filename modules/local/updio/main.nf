process UPDIO {
    tag "$meta.id"
    label 'process_single'

    container "cmgg/updio:1.0.0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(cnv)

    output:
    tuple val(meta), path("${prefix}"), emit: bed
    path  "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def common_cnv_file = cnv ? "--common_cnv_file $cnv" : "--common_cnv_file /usr/local/lib/updio/sample_data/common_dels_1percent_liftover.tsv"
    def VERSION = "1.0.0"

    """
    UPDio \\
        --multisample_vcf $vcf \\
        --output_path $prefix \\
        $common_cnv_file \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UPDio: $VERSION
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    def VERSION = "1.0.0"
    """
    mkdir $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        UPDio: $VERSION
    END_VERSIONS
    """
}
