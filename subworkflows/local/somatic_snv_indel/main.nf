//
// Somatic SNV/indel calling with DEEPSOMATIC, annotation with VEP, and optional somatic phasing
//

include { DEEPSOMATIC               } from '../../../modules/nf-core/deepsomatic/main'
include { ENSEMBLVEP_DOWNLOAD       } from '../../../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP            } from '../../../modules/nf-core/ensemblvep/vep/main'
include { HIPHASE_SOMATIC           } from '../../../modules/local/hiphase_somatic/main'
include { UNZIPFILES                } from '../../../modules/nf-core/unzipfiles/main'

workflow SOMATIC_SNV_INDEL {

    take:
    ch_tn_bam_pairs     // channel: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_genome_fasta     // channel: [ meta, fasta ]
    ch_genome_fai       // channel: [ meta, fai ]
    ch_vcf_tumors       // channel: [ meta, vcf, tbi ] - germline VCFs for tumor samples
    ch_bam_tumors       // channel: [ meta, bam, bai ] - tumor BAMs
    skip_vep            // boolean: skip VEP annotation
    skip_vep_download   // boolean: skip VEP cache download (use existing cache)
    skip_somatic_hiphase // boolean: skip somatic variant phasing

    main:
    ch_versions = Channel.empty()

    //
    // DEEPSOMATIC: Somatic SNV/indel calling
    //
    ch_deepsomatic_interval = channel.of( [[:], []] )
    ch_deepsomatic_gzi      = channel.of( [[:], []] )

    DEEPSOMATIC(
        ch_tn_bam_pairs,
        ch_deepsomatic_interval,
        ch_genome_fasta,
        ch_genome_fai,
        ch_deepsomatic_gzi
    )
    ch_versions = ch_versions.mix(DEEPSOMATIC.out.versions)

    ch_somatic_vcf     = DEEPSOMATIC.out.vcf
    ch_somatic_vcf_tbi = DEEPSOMATIC.out.vcf_tbi

    //
    // VEP: Variant annotation (optional)
    //
    ch_vep_vcf = Channel.empty()
    ch_vep_tab = Channel.empty()

    if (!skip_vep) {
        // Validate VEP parameters
        if (!params.vep_assembly || !params.vep_species || !params.vep_cache_version) {
            log.warn "VEP annotation requires vep_assembly, vep_species, and vep_cache_version parameters. Skipping VEP."
        } else {
            vep_cache_path = params.vep_cache ?: []

            // Download VEP cache if needed
            if (!skip_vep_download && !params.vep_cache) {
                ch_vep_download = Channel.of(
                    [[:], params.vep_assembly, params.vep_species, params.vep_cache_version]
                )

                ENSEMBLVEP_DOWNLOAD(ch_vep_download)
                ch_versions = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

                vep_cache_path = ENSEMBLVEP_DOWNLOAD.out.cache
                    .map { meta, path_prefix -> [path_prefix] }
            }

            // Prepare VEP input
            ch_vep_input = ch_somatic_vcf
                .map { meta, vcf_gz -> [ meta, vcf_gz, [] ] }

            ENSEMBLVEP_VEP(
                ch_vep_input,
                params.vep_assembly,
                params.vep_species,
                params.vep_cache_version,
                vep_cache_path,
                ch_genome_fasta,
                []
            )
            ch_vep_vcf = ENSEMBLVEP_VEP.out.vcf
            ch_vep_tab = ENSEMBLVEP_VEP.out.tab
            ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
        }
    }

    //
    // HIPHASE_SOMATIC: Somatic variant phasing (optional)
    //
    ch_somatic_phased_bam_bai = Channel.empty()
    ch_somatic_phased_vcf     = Channel.empty()

    if (!skip_somatic_hiphase) {
        // Prepare tumor germline VCF channel with patient_tumor_id as key
        ch_tumor_hiphase_vcf = ch_vcf_tumors
            .map { meta, vcf, tbi ->
                [meta.id, meta, vcf, tbi]
            }

        // Prepare somatic VCF channel with patient_tumor_id as key
        ch_somatic_hiphasing_vcf = ch_somatic_vcf
            .join(ch_somatic_vcf_tbi)
            .map { pair_meta, snv_vcf, vcf_tbi ->
                def patient_tumor_id = "${pair_meta.patient}_${pair_meta.tumor_id}"
                [ patient_tumor_id, pair_meta, snv_vcf, vcf_tbi]
            }

        // Prepare tumor BAM channel with patient_tumor_id as key
        ch_somatic_hiphasing_bam_bai = ch_bam_tumors
            .map { meta, bam, bai ->
                [meta.id, meta, bam, bai]
            }

        // Combine channels
        ch_somatic_hiphasing_combine = ch_tumor_hiphase_vcf
            .combine(ch_somatic_hiphasing_bam_bai, by: 0)
            .combine(ch_somatic_hiphasing_vcf, by: 0)
            .multiMap { meta_id, meta, vcf, tbi, meta2, bam, bai, meta3, somatic_vcf, somatic_tbi ->
                bam_bai:          [meta2, bam, bai]
                vcf_tbi:          [meta, vcf, tbi]
                somatic_vcf_tbi:  [meta3, somatic_vcf, somatic_tbi]
            }

        HIPHASE_SOMATIC(
            ch_somatic_hiphasing_combine.vcf_tbi,
            ch_somatic_hiphasing_combine.bam_bai,
            ch_genome_fasta,
            ch_somatic_hiphasing_combine.somatic_vcf_tbi
        )

        ch_somatic_phased_bam_bai = HIPHASE_SOMATIC.out.bam.join(HIPHASE_SOMATIC.out.bai)
        // ch_somatic_phased_vcf     = HIPHASE_SOMATIC.out.vcf
        ch_versions = ch_versions.mix(HIPHASE_SOMATIC.out.versions)
    }

    emit:
    vcf                = ch_somatic_vcf           // channel: [ pair_meta, vcf.gz ]
    vcf_tbi            = ch_somatic_vcf_tbi       // channel: [ pair_meta, vcf.gz.tbi ]
    vep_vcf            = ch_vep_vcf               // channel: [ pair_meta, annotated.vcf ]
    vep_tab            = ch_vep_tab               // channel: [ pair_meta, vep.tab ]
    phased_bam         = ch_somatic_phased_bam_bai // channel: [ meta, phased_bam, bai ]
    // phased_vcf         = ch_somatic_phased_vcf    // channel: [ meta, phased_vcf ]
    versions           = ch_versions              // channel: [ versions.yml ]
}
