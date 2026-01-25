process PURPLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-purple:4.0.2--hdfd78af_0' :
        'biocontainers/hmftools-purple:4.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(amber), path(cobalt), path(sv_hard_vcf), path(sv_hard_vcf_index), path(sv_soft_vcf), path(sv_soft_vcf_index), path(smlv_tumor_vcf), path(smlv_normal_vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    val(genome_ver)
    path(gc_profile)
    path(sage_known_hotspots_somatic)
    path(sage_known_hotspots_germline)
    path(driver_gene_panel)
    path(ensembl_data_dir)
    path(germline_del_freq)

    output:
    tuple val(meta), path("purple/"), emit: purple_dir
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_sample = meta.tumor_id ?: meta.id
    def normal_sample = meta.normal_id ?: "${meta.id}_normal"

    def smlv_tumor_arg = smlv_tumor_vcf ? "-somatic_vcf prepared_somatic.vcf.gz" : ""
    def smlv_normal_arg = smlv_normal_vcf ? "-germline_vcf prepared_germline.vcf.gz" : ""
    def germline_del_freq_arg = germline_del_freq ? "-germline_del_freq_file ${germline_del_freq}" : ""

    """
    # For provided smlv VCFs, filter records that do not contain the required FORMAT/AD field
    if [ -n "${smlv_tumor_vcf}" ]; then
        bcftools filter -Oz -e 'FORMAT/AD[*]="."' "${smlv_tumor_vcf}" > prepared_somatic.vcf.gz
    fi

    if [ -n "${smlv_normal_vcf}" ]; then
        bcftools filter -Oz -e 'FORMAT/AD[*]="."' "${smlv_normal_vcf}" > prepared_germline.vcf.gz
    fi

    # Run PURPLE  # -reference        -sv_recovery_vcf "${sv_soft_vcf}" \\
    #        -structural_vcf "${sv_hard_vcf}" \\   -gc_profile     -run_drivers \\

    purple \\
        $args \\
        -tumor "${tumor_sample}" \\
        -reference "${normal_sample}" \\
        $smlv_tumor_arg \\
        $smlv_normal_arg \\
        -amber "${amber}" \\
        -cobalt "${cobalt}" \\
        -output_dir purple/ \\
        -gc_profile "${gc_profile}" \\
        -driver_gene_panel "${driver_gene_panel}" \\
        -ensembl_data_dir "${ensembl_data_dir}" \\
        -somatic_hotspots "${sage_known_hotspots_somatic}" \\
        -germline_hotspots "${sage_known_hotspots_germline}" \\
        $germline_del_freq_arg \\
        -ref_genome "${fasta}" \\
        -ref_genome_version "${genome_ver}" \\
        -threads ${task.cpus}

    # PURPLE can fail silently, check that at least the PURPLE SV VCF is created
    # if [ ! -s "purple/${tumor_sample}.purple.sv.vcf.gz" ]; then
    #    echo "ERROR: PURPLE output file not created or empty" >&2
    #    exit 1
    # fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: 4.0.2
        bcftools: 1.19
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor_sample = meta.tumor_id ?: meta.id
    """
    mkdir -p purple/
    touch purple/${tumor_sample}.purple.sv.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purple: 4.0.2
        bcftools: 1.19
    END_VERSIONS
    """
}
