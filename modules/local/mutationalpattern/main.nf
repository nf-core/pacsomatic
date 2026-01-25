process MUTATIONALPATTERN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-mutationalpatterns_r-biocmanager:81c7d2739eed68a7':
        'community.wave.seqera.io/library/bioconductor-mutationalpatterns_r-biocmanager:b3aab5922ca450c5' }"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), val(ref_genome)
    val(max_delta)

    output:
    tuple val(meta), path("*.mut_sigs.tsv")             , emit: mut_sig
    tuple val(meta), path("*.reconstructed_sigs.tsv")   , emit: recon_sig
    tuple val(meta), path("*.type_occurences.tsv")      , emit: type
    tuple val(meta), path("*.mut_sigs_bootstrapped.tsv"), emit: mut_sig_boot
    tuple val(meta), path("*.mutation_profile.pdf")     , emit: mut_profile
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'mutational_pattern.R'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mut_sigs.tsv
    touch ${prefix}.reconstructed_sigs.tsv
    touch ${prefix}.type_occurences.tsv
    touch ${prefix}.mut_sigs_bootstrapped.tsv
    touch ${prefix}.mutation_profile.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r: \$(R --version | grep -oP "\\d+\\.\\d+\\.\\d+")
        mutationalpattern: \$(Rscript -e "library(MutationalPatterns); cat(as.character(packageVersion('MutationalPatterns')), '\\n', sep='')")
        cowplot: \$(Rscript -e "library(cowplot); cat(as.character(packageVersion('cowplot')), '\\n', sep='')")
        ggplot2: \$(Rscript -e "library(ggplot2); cat(as.character(packageVersion('ggplot2')), '\\n', sep='')")
    END_VERSIONS
    """
}
