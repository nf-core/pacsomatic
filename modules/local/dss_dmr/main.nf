process DSS_DMR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-bsseq_bioconductor-dss_r-data.table:79ccb00e9e587149':
        'community.wave.seqera.io/library/bioconductor-bsseq_bioconductor-dss_r-data.table:49f7b3fac7d51eed' }"

    input:
    tuple val(meta), path(tumor_bed), path(normal_bed)

    output:
    tuple val(meta), path("*.dmr.tsv")  , emit: dmr
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'dss_dmr.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.dmr.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | grep -oP "\\d+\\.\\d+\\.\\d+")
        dss: \$(Rscript -e "library(DSS); cat(as.character(packageVersion('DSS')), '\\n', sep='')")
        bsseq: \$(Rscript -e "library(bsseq); cat(as.character(packageVersion('bsseq')), '\\n', sep='')")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')), '\\n', sep='')")
    END_VERSIONS
    """
}
