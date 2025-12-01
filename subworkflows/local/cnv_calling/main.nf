//
// Copy number variant calling with CNVKIT
//

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL  } from '../../../modules/nf-core/cnvkit/call/main'

workflow CNV_CALLING {

    take:
    ch_tn_bam_pairs     // channel: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_genome_fasta     // channel: [ meta, fasta ]
    ch_genome_fai       // channel: [ meta, fai ]

    main:
    ch_versions = Channel.empty()

    //
    // CNVKIT_BATCH: Initial CNV calling
    //
    // Reorder to [meta, tumor_bam, normal_bam]
    ch_cnvkit_input = ch_tn_bam_pairs
        .map { pair_meta, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ->
            [pair_meta, tumor_bam, normal_bam]
        }

    // Prepare optional inputs
    ch_cnvkit_targets = params.cnvkit_targets ?
        channel.of([[:], file(params.cnvkit_targets, checkIfExists: true)]) :
        channel.of([[:], []])

    ch_cnvkit_reference = params.cnvkit_reference ?
        channel.of([[:], file(params.cnvkit_reference, checkIfExists: true)]) :
        channel.of([[:], []])

    CNVKIT_BATCH(
        ch_cnvkit_input,
        ch_genome_fasta,
        ch_genome_fai,
        ch_cnvkit_targets,
        ch_cnvkit_reference,
        [:]
    )
    ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

    //
    // CNVKIT_CALL: Refine CNV calls
    //
    // Optional: Merge germline VCFs for refinement
    ch_merged_germline_vcf = params.cnvkit_vcf ?
        channel.of(file(params.cnvkit_vcf, checkIfExists: true)) :
        channel.of([[]])

    // Extract call_cns from the CNS output (index 1 of the list)
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
