process BCFTOOLS_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.${extension}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def last_args = args3 ?: args2 ?: args
    def prefix = task.ext.prefix ?: "${meta.id}"

    extension = last_args.contains("--output-type b") || last_args.contains("-Ob")
        ? "bcf.gz"
        : last_args.contains("--output-type u") || last_args.contains("-Ou")
            ? "bcf"
            : last_args.contains("--output-type z") || last_args.contains("-Oz")
                ? "vcf.gz"
                : last_args.contains("--output-type v") || last_args.contains("-Ov")
                    ? "vcf"
                    : "vcf"

    if ("${vcf}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    def filter_2 = args2 ? "| bcftools filter --threads ${task.cpus} ${args2}" : ""
    def filter_3 = args3 ? "| bcftools filter --threads ${task.cpus} ${args3}" : ""

    """
    bcftools filter \\
        --output ${prefix}.${extension} \\
        --threads ${task.cpus} \\
        ${args} \\
        ${vcf} \\
        ${filter_2} \\
        ${filter_3} \\
        --output ${prefix}.${extension}
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def last_args = args3 ?: args2 ?: args
    def prefix = task.ext.prefix ?: "${meta.id}"
    extension = last_args.contains("--output-type b") || last_args.contains("-Ob")
        ? "bcf.gz"
        : last_args.contains("--output-type u") || last_args.contains("-Ou")
            ? "bcf"
            : last_args.contains("--output-type z") || last_args.contains("-Oz")
                ? "vcf.gz"
                : last_args.contains("--output-type v") || last_args.contains("-Ov")
                    ? "vcf"
                    : "vcf"
    def index = last_args.contains("--write-index=tbi") || last_args.contains("-W=tbi")
        ? "tbi"
        : last_args.contains("--write-index=csi") || last_args.contains("-W=csi")
            ? "csi"
            : last_args.contains("--write-index") || last_args.contains("-W")
                ? "csi"
                : ""
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def create_index = extension.endsWith(".gz") && index.matches("csi|tbi") ? "touch ${prefix}.${extension}.${index}" : ""

    if ("${vcf}" == "${prefix}.${extension}") {
        error("Input and output names are the same, set prefix in module configuration to disambiguate!")
    }

    """
    ${create_cmd} ${prefix}.${extension}
    ${create_index}
    """
}
