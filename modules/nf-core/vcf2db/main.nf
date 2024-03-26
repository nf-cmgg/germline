process VCF2DB {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "cmgg/vcf2db:2020.02.24"

    input:
    tuple val(meta), path(vcf), path(ped)

    output:
    tuple val(meta), path("*.db") , emit: db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2020.02.24" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    vcf2db.py \\
        $vcf \\
        $ped \\
        ${prefix}.db \\
        $args

    sqlite3 ${prefix}.db 'CREATE INDEX idx_variant_impacts_id ON variant_impacts (variant_id)' && \\
    sqlite3 ${prefix}.db 'ALTER TABLE variants ADD COLUMN tags varchar(255)' && \\
    sqlite3 ${prefix}.db 'ALTER TABLE variants ADD COLUMN tags_user varchar(255)' && \\
    sqlite3 ${prefix}.db 'ALTER TABLE variants ADD COLUMN notes varchar(255)' && \\
    sqlite3 ${prefix}.db 'ALTER TABLE variants ADD COLUMN notes_user varchar(255)'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2db: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2020.02.24" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2db: $VERSION
    END_VERSIONS
    """
}
