process SVPACK_MATCH {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_pysam:0702f875be118f5c':
        'community.wave.seqera.io/library/htslib_pysam:ba875304e3f7d749' }"

    input:
    tuple val(meta), path(vcf_a)
    tuple val(meta2), path(vcf_b)

    output:
    tuple val(meta), path("*.matched.vcf"), emit: vcf
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2025.10.31" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    svpack match \\
        $vcf_a \\
        $vcf_b \\
        $args \\
        > ${prefix}.matched.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svpack: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2025.10.31" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.matched.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svpack: $VERSION
    END_VERSIONS
    """
}
