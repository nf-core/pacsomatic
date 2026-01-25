//
// Tumor clonality analysis using Hartwig Medical Foundation tools: AMBER, COBALT, and PURPLE
//

include { AMBER  } from '../../../modules/local/amber/main'
include { COBALT } from '../../../modules/local/cobalt/run/main'
include { PURPLE } from '../../../modules/local/purple/main'

workflow TUMOR_CLONALITY {

    take:
    ch_tn_bam_pairs                // channel: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_genome_fasta                // channel: [ meta, fasta ]
    ch_genome_fai                  // channel: [ meta, fai ]
    ch_heterozygous_sites          // channel: [ vcf(.gz) ]
    ch_target_regions_bed          // channel: [ bed ]
    ch_gc_profile                  // channel: [ cnp ]
    ch_diploid_regions             // channel: [ bed(.gz) ]
    ch_target_region_normalisation // channel: [ tsv ]
    ch_known_hotspots_somatic      // channel: [ vcf(.gz) ]
    ch_known_hotspots_gemrline     // channel: [ vcf(.gz) ]
    ch_driver_gene_panel           // channel: [ tsv ]
    ch_ensembl_data_dir            // channel: [ dir ]

    main:
    ch_versions = Channel.empty()

    // Validate required parameters
    def required_params = [
        'heterozygous_sites',
        'gc_profile',
        'driver_gene_panel',
        'known_hotspots_somatic',
        'known_hotspots_germline',
        'ensembl_data_dir'
    ]

    def missing_params = required_params.findAll { !params[it] }
    if (missing_params) {
        log.error "Tumor clonality analysis requires the following parameters: ${missing_params.join(', ')}"
        log.error "Please provide these parameters or use --skip_tumor_clonality to skip this analysis."
        exit 1
    }

    //
    // AMBER: Allelic ratio estimation
    //
    ch_amber_input = ch_tn_bam_pairs
        .map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
            // [meta, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai]
            [meta, tumor_bam, normal_bam, [], tumor_bai, normal_bai, []]
        }

    // Prepare optional parameters

    AMBER(
        ch_amber_input,
        'V38',
        ch_heterozygous_sites,
        ch_target_regions_bed,
        []
    )
    ch_versions = ch_versions.mix(AMBER.out.versions)

    // Add pair ID as key for joining
    ch_amber_dir = AMBER.out.amber_dir
        .map { pair_meta, amber_dir ->
            [pair_meta.id, pair_meta, amber_dir]
        }

    //
    // COBALT: Read depth ratio estimation
    //
    ch_cobalt_input = ch_tn_bam_pairs
        .map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
            [meta, tumor_bam, normal_bam, tumor_bai, normal_bai]
        }

    // Prepare optional parameters

    COBALT(
        ch_cobalt_input,
        ch_gc_profile,
        ch_diploid_regions,
        ch_target_region_normalisation,
        [:]
    )
    ch_versions = ch_versions.mix(COBALT.out.versions)

    // Add pair ID as key for joining
    ch_cobalt_dir = COBALT.out.cobalt_dir
        .map { pair_meta, cobalt_dir ->
            [pair_meta.id, pair_meta, cobalt_dir]
        }

    //
    // PURPLE: Purity and ploidy estimation
    //
    // Join AMBER and COBALT outputs
    ch_purple_input = AMBER.out.amber_dir
        .join(COBALT.out.cobalt_dir)
        .map { meta, amber_dir, cobalt_dir ->
            // [meta, amber_dir, cobalt_dir, sv_hard_vcf, sv_hard_vcf_index,
            //  sv_soft_vcf, sv_soft_vcf_index, smlv_tumor_vcf, smlv_normal_vcf]
            [meta, amber_dir, cobalt_dir, [], [], [], [], [], []]
        }


    PURPLE(
        ch_purple_input,
        ch_genome_fasta,
        ch_genome_fai,
        [[:], []],
        '38',
        ch_gc_profile,
        ch_known_hotspots_somatic,
        ch_known_hotspots_gemrline,
        ch_driver_gene_panel,
        ch_ensembl_data_dir,
        []
    )
    ch_versions = ch_versions.mix(PURPLE.out.versions)

    emit:
    amber_dir   = AMBER.out.amber_dir       // channel: [ pair_meta, amber_dir ]
    cobalt_dir  = COBALT.out.cobalt_dir     // channel: [ pair_meta, cobalt_dir ]
    purple_dir  = PURPLE.out.purple_dir     // channel: [ pair_meta, purple_dir ]
    // purity_tsv  = PURPLE.out.purity_tsv     // channel: [ pair_meta, purity.tsv ]
    versions    = ch_versions               // channel: [ versions.yml ]
}
