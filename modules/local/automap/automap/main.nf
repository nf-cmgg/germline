process AUTOMAP_AUTOMAP {
    tag "$meta.id"
    label 'process_single'

    container "cmgg/automap:1.0.0"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple val(meta2), path(repeats)
    tuple val(meta3), path(panel)
    val(genome)

    output:
    tuple val(meta), path("${prefix}"), emit: automap
    path  "versions.yml"              , emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def panel_file = panel ? "--panel $panel" : "--panel /usr/local/lib/automap/Resources/Biomodule_20220808_all_genes_hg38.txt"
    def hg_genome = genome ?: "hg38"
    def VERSION = "1.0.0"

    """
    automap \\
        --vcf $vcf \\
        --genome $hg_genome \\
        --out $prefix/ \\
        --repeats $repeats \\
        $panel_file \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        automap: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def panel_name = args.contains("--panelname") ? args.split("--panelname")[-1].trim().split(" ")[0] : ""
    prefix = task.ext.prefix ?: "${meta.id}"

    def create_outputs = meta.family_count > 1 ? (1..meta.family_count).collect { number ->
        def cmd_prefix = "touch ${prefix}/sample${number}"
        [
            "mkdir ${prefix}/sample${number}",
            "${cmd_prefix}/sample${number}.HomRegions.pdf",
            "${cmd_prefix}/sample${number}.HomRegions.${panel_name}.tsv",
            "${cmd_prefix}/sample${number}.HomRegions.tsv",
            "${cmd_prefix}/sample${number}.HomRegions.strict.${panel_name}.tsv"
        ].join(" && ")
    }.join(" && ") : [
        "touch ${prefix}/${meta.id}.HomRegions.pdf",
        "touch ${prefix}/${meta.id}.HomRegions.${panel_name}.tsv",
        "touch ${prefix}/${meta.id}.HomRegions.tsv",
        "touch ${prefix}/${meta.id}.HomRegions.strict.${panel_name}.tsv"
    ].join(" && ")

    def VERSION = "1.0.0"
    """
    mkdir $prefix
    ${create_outputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        automap: $VERSION
    END_VERSIONS
    """
}
