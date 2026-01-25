process SVPACK_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/htslib_pysam:0702f875be118f5c':
        'community.wave.seqera.io/library/htslib_pysam:ba875304e3f7d749' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: vcf
    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2025.10.31" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    svpack filter \\
        $vcf \\
        $args \\
        > ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svpack: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "2025.10.31" // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        svpack: $VERSION
    END_VERSIONS
    """
}
