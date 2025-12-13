//
// Somatic structural variant calling with SEVERUS, filtering with SVPACK, and annotation with ANNOTSV
//

include { SEVERUS                    } from '../../../modules/nf-core/severus/main'
include { SVPACK_ANNOTATE            } from '../svpack_annotate/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../../../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../../../modules/nf-core/annotsv/annotsv/main'
include { TABIX_BGZIPTABIX as TABIX_SV_VCF } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow SOMATIC_SV {

    take:
    ch_tn_bam_pairs        // channel: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_severus_trf_bed     // channel: [ meta, bed ]
    ch_svpack_ctrl_vcf     // channel: [ meta, vcf(.gz) ]
    ch_svpack_ref_gff      // channel: [ meta, gff(.gz) ]
    ch_annotsv_cache       // channel: [ meta, annotsv_cache ]
    skip_svpack            // boolean: skip SVPACK filtering
    skip_annotsv           // boolean: skip ANNOTSV annotation
    skip_annotsv_install   // boolean: skip ANNOTSV cache installation (use existing cache)

    main:
    ch_versions = Channel.empty()

    //
    // SEVERUS: Somatic SV calling
    //

    // Reorder BAMs: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    ch_severus_input = ch_tn_bam_pairs
        .map { pair_meta, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ->
            [pair_meta, tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai, []]
        }

    SEVERUS(ch_severus_input, ch_severus_trf_bed)
    ch_versions = ch_versions.mix(SEVERUS.out.versions)

    ch_sv_vcf = SEVERUS.out.somatic_vcf

    //
    // SVPACK: SV filtering and annotation (optional)
    //
    ch_sv_filtered_vcf = ch_sv_vcf

    if (!skip_svpack) {
        // Validate SVPACK parameters
        if (!params.svpack_ctrl_vcf || !params.svpack_ref_gff) {
            log.warn "SVPACK filtering requires svpack_control_vcf and svpack_ref_gff parameters. Skipping SVPACK."
        } else {
            SVPACK_ANNOTATE(ch_sv_vcf, ch_svpack_ctrl_vcf, ch_svpack_ref_gff)
            ch_sv_filtered_vcf = SVPACK_ANNOTATE.out.tagged_vcf
            ch_versions = ch_versions.mix(SVPACK_ANNOTATE.out.versions)
        }
    }

    // Compress and index SV VCF
    TABIX_SV_VCF(ch_sv_filtered_vcf)
    ch_versions = ch_versions.mix(TABIX_SV_VCF.out.versions)

    //
    // ANNOTSV: SV annotation (optional)
    //
    ch_annotsv_tsv = Channel.empty()

    if (!skip_annotsv) {
        // Validate ANNOTSV parameters
        if (!params.annotsv_cache && skip_annotsv_install) {
            log.warn "ANNOTSV annotation requires annotsv_cache parameter or skip_annotsv_install=false. Skipping ANNOTSV."
        } else {
            // Install ANNOTSV annotations if needed
            if (!skip_annotsv_install) {
                ANNOTSV_INSTALLANNOTATIONS()
                ch_versions = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
                ch_annotsv_cache = ANNOTSV_INSTALLANNOTATIONS.out.annotations
                    .map { AnnotSV_annotations -> [[:], AnnotSV_annotations] }
            }

            // Prepare input for ANNOTSV
            ch_annotsv_input = TABIX_SV_VCF.out.gz_tbi
                .map { meta, vcf, tbi ->
                    [meta, vcf, tbi, []]  // [meta, sv_vcf, sv_vcf_tbi, candidate_small_variants]
                }

            ANNOTSV_ANNOTSV(
                ch_annotsv_input,
                ch_annotsv_cache,
                [[:], []],
                [[:], []],
                [[:], []]
            )
            ch_annotsv_tsv = ANNOTSV_ANNOTSV.out.tsv
            ch_versions = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions)
        }
    }

    emit:
    vcf         = ch_sv_vcf         // channel: [ pair_meta, vcf ]
    annotsv_tsv = ch_annotsv_tsv    // channel: [ pair_meta, tsv ]
    versions    = ch_versions       // channel: [ versions.yml ]
}
