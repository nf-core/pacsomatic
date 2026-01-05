//
// Methylation analysis with PBCPGTOOLS and differential methylation region (DMR) detection
//

include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES as PBCPGTOOLS_NORMAL } from '../../../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES as PBCPGTOOLS_TUMOR  } from '../../../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { DSS_DMR                                                } from '../../../modules/local/dss_dmr/main'
include { ANNOTATR_DMR                                           } from '../../../modules/local/annotatr_dmr/main'

workflow METHYLATION_ANALYSIS {

    take:
    ch_phased_normal_bam_bai  // channel: [ meta, bam, bai ] - phased normal samples
    ch_phased_tumor_bam_bai   // channel: [ meta, bam, bai ] - phased tumor samples
    skip_dmr                  // boolean: skip DMR detection
    skip_dmr_anno             // boolean: skip DMR annotation

    main:
    ch_versions = Channel.empty()

    //
    // PBCPGTOOLS: CpG methylation calling on normal samples
    //
    PBCPGTOOLS_NORMAL(ch_phased_normal_bam_bai)
    ch_versions = ch_versions.mix(PBCPGTOOLS_NORMAL.out.versions.first())

    // Add patient as key for pairing
    ch_normal_cpg_bed = PBCPGTOOLS_NORMAL.out.combined_bed
        .map { meta, cpg_bed ->
            [meta.id, meta, cpg_bed]
        }

    //
    // PBCPGTOOLS: CpG methylation calling on tumor samples
    //
    PBCPGTOOLS_TUMOR(ch_phased_tumor_bam_bai)
    ch_versions = ch_versions.mix(PBCPGTOOLS_TUMOR.out.versions.first())

    // Add patient as key for pairing
    ch_tumor_cpg_bed = PBCPGTOOLS_TUMOR.out.combined_bed
        .map { pair_meta, cpg_bed ->
            def patient_normal_id = "${pair_meta.patient}_${pair_meta.normal_id}"
            [patient_normal_id, pair_meta, cpg_bed]
        }

    //
    // DMR: Differential methylation region detection (optional)
    //
    ch_dmr_tsv = Channel.empty()
    ch_dmr_annotated = Channel.empty()

    if (!skip_dmr) {
        // Pair tumor and normal CpG BED files by patient
        ch_tn_pair_dmr = ch_tumor_cpg_bed
            .combine(ch_normal_cpg_bed, by: [0])
            .map { patient_normal_id, pair_meta, tumor_bed, normal_meta, normal_bed ->
                [pair_meta, tumor_bed, normal_bed]
            }

        DSS_DMR(ch_tn_pair_dmr)
        ch_dmr_tsv = DSS_DMR.out.dmr
        ch_versions = ch_versions.mix(DSS_DMR.out.versions.first())

        //
        // ANNOTATR: DMR annotation (optional)
        //
        if (!skip_dmr_anno) {
            ANNOTATR_DMR(ch_dmr_tsv)
            ch_dmr_annotated = ANNOTATR_DMR.out.summary
            ch_versions = ch_versions.mix(ANNOTATR_DMR.out.versions.first())
        }
    }

    emit:
    normal_cpg      = PBCPGTOOLS_NORMAL.out.combined_bed  // channel: [ meta, bed ]
    tumor_cpg       = PBCPGTOOLS_TUMOR.out.combined_bed   // channel: [ meta, bed ]
    dmr             = ch_dmr_tsv                          // channel: [ meta, tsv ]
    dmr_annotated   = ch_dmr_annotated                    // channel: [ meta, annotated ]
    versions        = ch_versions                         // channel: [ versions.yml ]
}
