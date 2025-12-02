//
// Copy number variant calling with CNVKIT
//

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL  } from '../../../modules/nf-core/cnvkit/call/main'

workflow CNV_CALLING {

    take:
    ch_tn_bam_pairs        // channel: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_genome_fasta        // channel: [ meta, fasta ]
    ch_genome_fai          // channel: [ meta, fai ]
    ch_cnv_target_bed      // channel: [ meta, bed ]
    ch_cnv_reference       // channel: [ meta, cnn ]
    ch_merged_germline_vcf // channel: [ meta, vcf ]

    main:
    ch_versions = Channel.empty()

    //
    // CNVKIT_BATCH: Initial CNV calling
    //
    // Reorder to [meta, tumor_bam, normal_bam]
    ch_cnvkit_input = ch_tn_bam_pairs
        .map { pair_meta, normal_bam, _normal_bam_bai, tumor_bam, _tumor_bam_bai ->
            [pair_meta, tumor_bam, normal_bam]
        }

    CNVKIT_BATCH(
        ch_cnvkit_input,
        ch_genome_fasta,
        ch_genome_fai,
        ch_cnv_target_bed,
        ch_cnv_reference,
        false
    )
    ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())

    //
    // CNVKIT_CALL: Refine CNV calls
    //

    // Extract call_cns from the CNS output (index 1 of the list)
    // Optional: Merge germline VCFs for refinement
    ch_cnvkit_call_input = CNVKIT_BATCH.out.cns
        .combine(ch_merged_germline_vcf)
        .map { cnvkit_meta, cns_list, merge_vcf ->
            [cnvkit_meta, cns_list[1], merge_vcf]  // [meta, call_cns, merge_vcf]
        }

    CNVKIT_CALL(ch_cnvkit_call_input)
    ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

    emit:
    cns         = CNVKIT_BATCH.out.cns      // channel: [ meta, [bintest_cns, call_cns, cns] ]
    cnr         = CNVKIT_BATCH.out.cnr      // channel: [ meta, cnr ]
    called_cns  = CNVKIT_CALL.out.cns       // channel: [ meta, called.cns ]
    versions    = ch_versions               // channel: [ versions.yml ]
}
