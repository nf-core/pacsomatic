/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SVPACK_FILTER      } from '../../../modules/local/svpack/filter/main'
include { SVPACK_MATCH       } from '../../../modules/local/svpack/match/main'
include { SVPACK_CONSEQUENCE } from '../../../modules/local/svpack/consequence/main'
include { SVPACK_TAGZYGOSITY } from '../../../modules/local/svpack/tagzygosity/main'
include { TABIX_BGZIPTABIX   } from '../../../modules/nf-core/tabix/bgziptabix/main'

workflow SVPACK_ANNOTATE {

    take:
    ch_vcf          // channel: [ val(meta), path(vcf) ]
    ch_vcf_match    // channel: [ val(meta2), path(vcf_b) ] - reference VCF for matching
    ch_genes_gff    // channel: [ val(meta3), path(gff) ] - genes annotation GFF file

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Filter structural variants
    //
    SVPACK_FILTER (
        ch_vcf
    )
    ch_versions = ch_versions.mix(SVPACK_FILTER.out.versions.first())

    //
    // MODULE: Match against reference VCF
    //
    SVPACK_MATCH (
        SVPACK_FILTER.out.vcf,
        ch_vcf_match
    )
    ch_versions = ch_versions.mix(SVPACK_MATCH.out.versions.first())

    //
    // MODULE: Annotate consequences on genes
    //
    SVPACK_CONSEQUENCE (
        SVPACK_MATCH.out.vcf,
        ch_genes_gff
    )
    ch_versions = ch_versions.mix(SVPACK_CONSEQUENCE.out.versions.first())

    //
    // MODULE: Tag zygosity information
    //
    SVPACK_TAGZYGOSITY (
        SVPACK_CONSEQUENCE.out.vcf
    )
    ch_versions = ch_versions.mix(SVPACK_TAGZYGOSITY.out.versions.first())

    //
    // MODULE: Compress and index final VCF
    //
    TABIX_BGZIPTABIX (
        SVPACK_TAGZYGOSITY.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    emit:
    filtered_vcf    = SVPACK_FILTER.out.vcf        // channel: [ val(meta), path(vcf) ]
    matched_vcf     = SVPACK_MATCH.out.vcf         // channel: [ val(meta), path(vcf) ]
    consequence_vcf = SVPACK_CONSEQUENCE.out.vcf   // channel: [ val(meta), path(vcf) ]
    tagged_vcf      = SVPACK_TAGZYGOSITY.out.vcf   // channel: [ val(meta), path(vcf) ]
    vcf_gz          = TABIX_BGZIPTABIX.out.gz_tbi  // channel: [ val(meta), path(vcf.gz), path(vcf.gz.tbi) ]
    versions        = ch_versions                  // channel: [ path(versions.yml) ]
}
