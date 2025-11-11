process ANNOTATR_DMR {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-annotatr_bioconductor-genomicranges_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-data.table:30511d9a64426ffb':
        'community.wave.seqera.io/library/bioconductor-annotatr_bioconductor-genomicranges_bioconductor-org.hs.eg.db_bioconductor-txdb.hsapiens.ucsc.hg38.knowngene_r-data.table:e048f72fc22326e9' }"

    input:
    tuple val(meta), path(dmr)

    output:
    tuple val(meta), path("*_dmr_annotation_summary.tsv.gz"), emit: summary
    tuple val(meta), path("*_*_dmrs.tsv.gz")                , emit: annotated, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    template 'annotatr_dmr.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_dmr_annotation_summary.tsv.gz
    touch ${prefix}_hg38_genes_promoters_dmrs.tsv.gz
    touch ${prefix}_hg38_genes_1to5kb_dmrs.tsv.gz
    touch ${prefix}_hg38_genes_5UTRs_dmrs.tsv.gz
    touch ${prefix}_hg38_genes_exons_dmrs.tsv.gz
    touch ${prefix}_hg38_genes_introns_dmrs.tsv.gz
    touch ${prefix}_hg38_genes_3UTRs_dmrs.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | grep -oP "\\d+\\.\\d+\\.\\d+")
        annotatr: \$(Rscript -e "library(annotatr); cat(as.character(packageVersion('annotatr')), '\\n', sep='')")
        genomicranges: \$(Rscript -e "library(GenomicRanges); cat(as.character(packageVersion('GenomicRanges')), '\\n', sep='')")
        data.table: \$(Rscript -e "library(data.table); cat(as.character(packageVersion('data.table')), '\\n', sep='')")
    END_VERSIONS
    """
}
