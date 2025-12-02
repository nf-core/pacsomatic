//
// Germline variant calling with CLAIR3 and optional phasing with HIPHASE
//

include { CLAIR3  } from '../../../modules/nf-core/clair3/main'
include { HIPHASE } from '../../../modules/nf-core/hiphase/main'

workflow GERMLINE_CALLING_PHASING {

    take:
    ch_bam_bai           // channel: [ meta, bam, bai ]
    ch_genome_fasta      // channel: [ meta, fasta ]
    ch_genome_fai        // channel: [ meta, fai ]
    skip_hiphase         // boolean: skip HiPhase phasing

    main:
    ch_versions = Channel.empty()

    //
    // Validate CLAIR3 model configuration
    //
    if (params.clair3_model && params.clair3_model_path) {
        log.error "Two models specified ${params.clair3_model} and ${params.clair3_model_path}, specify one of them."
        exit 1
    }
    if (!params.clair3_model && !params.clair3_model_path) {
        log.error "No clair3 model is specified, use option params.clair3_model or params.clair3_model_path."
        exit 1
    }

    //
    // CLAIR3: Germline variant calling
    //
    ch_clair3_input = ch_bam_bai.map { meta, bam, bai ->
        def clair3_model_path = params.clair3_model_path ?
            file(params.clair3_model_path, checkIfExists:true, dir:true) : []
        [ meta, bam, bai, params.clair3_model, clair3_model_path, 'hifi' ]
    }

    CLAIR3 (
        ch_clair3_input,
        ch_genome_fasta,
        ch_genome_fai
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    // Join VCF and TBI
    ch_vcf_tbi = CLAIR3.out.vcf.join(CLAIR3.out.tbi)

    //
    // HIPHASE: Germline phasing (optional, only for normal samples)
    //
    ch_phased_bam_bai = Channel.empty()
    ch_phased_vcf_tbi = Channel.empty()

    if (!skip_hiphase) {
        // Filter for normal samples only (status == 0)
        ch_hiphase_vcf = ch_vcf_tbi
            .filter { meta, _vcf, _tbi -> meta.status == 0 }
            .map { meta, vcf, tbi ->
                [meta.id, meta, vcf, tbi]
            }

        ch_hiphase_bam_bai = ch_bam_bai
            .filter { meta, _bam, _bai -> meta.status == 0 }
            .map { meta, bam, bai ->
                [meta.id, meta, bam, bai]
            }

        // Combine VCF and BAM by sample ID
        ch_hiphase_prep = ch_hiphase_vcf
            .combine(ch_hiphase_bam_bai, by: [0])
            .multiMap { meta_id, meta, vcf, tbi, meta2, bam, bai ->
                vcf_tbi: [meta, vcf, tbi]
                bam_bai: [meta2, bam, bai]
            }

        HIPHASE (
            ch_hiphase_prep.vcf_tbi,
            ch_hiphase_prep.bam_bai,
            ch_genome_fasta
        )

        ch_phased_bam_bai = HIPHASE.out.bam.join(HIPHASE.out.bai)
        ch_phased_vcf_tbi = HIPHASE.out.vcf.join(HIPHASE.out.tbi)
        ch_versions = ch_versions.mix(HIPHASE.out.versions.first())
    }

    emit:
    vcf         = ch_vcf_tbi                // channel: [ meta, vcf, tbi ]
    phased_bam  = ch_phased_bam_bai         // channel: [ meta, bam, bai ] - phased normal samples
    phased_vcf  = ch_phased_vcf_tbi         // channel: [ meta, vcf, tbi ] - phased VCFs
    versions    = ch_versions               // channel: [ versions.yml ]
}
