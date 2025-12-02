//
// Mutational signature analysis with CHORD (HRD detection) and MUTATIONALPATTERN
//

include { CHORD               } from '../../../modules/local/chord/main'
include { MUTATIONALPATTERN   } from '../../../modules/local/mutationalpattern/main'
include { SAMTOOLS_DICT       } from '../../../modules/nf-core/samtools/dict/main'
include { UNZIPFILES          } from '../../../modules/nf-core/unzipfiles/main'

workflow SIGNATURE_ANALYSIS {

    take:
    ch_somatic_snv_vcf_gz   // channel: [ pair_meta, vcf.gz ]
    ch_somatic_sv_vcf       // channel: [ pair_meta, sv_vcf ]
    ch_genome_fasta         // channel: [ meta, fasta ]
    ch_genome_fai           // channel: [ meta, fai ]
    skip_chord              // boolean: skip CHORD HRD detection
    skip_mutationalpattern  // boolean: skip mutational pattern analysis

    main:
    ch_versions = Channel.empty()

    //
    // CHORD: Homologous recombination deficiency (HRD) detection
    //
    ch_chord_prediction = Channel.empty()

    if (!skip_chord) {
        // Unzip SNV VCF
        UNZIPFILES(ch_somatic_snv_vcf_gz)
        ch_versions = ch_versions.mix(UNZIPFILES.out.versions.first())

        // Prepare SNV VCF channel with pair_id as key
        ch_chord_somatic_snv_vcf = UNZIPFILES.out.files
            .map { pair_meta, snv_vcf ->
                [pair_meta.id, pair_meta, snv_vcf]
            }

        // Prepare SV VCF channel with pair_id as key
        ch_chord_somatic_sv_vcf = ch_somatic_sv_vcf
            .map { pair_meta, sv_vcf ->
                [pair_meta.id, pair_meta, sv_vcf]
            }

        // Combine SNV and SV VCFs
        ch_chord_input = ch_chord_somatic_snv_vcf
            .combine(ch_chord_somatic_sv_vcf, by: 0)
            .map { pair_id, pair_meta, snv_vcf, pair_meta2, sv_vcf ->
                def chord_meta = [
                    patient:    pair_meta.patient,
                    tumor_id:   pair_meta.tumor_id,
                    normal_id:  pair_meta.normal_id,
                    id:         pair_meta.id,
                    sample_id:  pair_meta.id
                ]
                [chord_meta, snv_vcf, sv_vcf]
            }

        // Prepare genome files
        chord_genome_fasta = ch_genome_fasta.map { meta, fasta -> [fasta] }
        chord_genome_fai   = ch_genome_fai.map { meta, fai -> [fai] }

        // Create genome dictionary
        SAMTOOLS_DICT(ch_genome_fasta)
        ch_versions = ch_versions.mix(SAMTOOLS_DICT.out.versions)

        chord_genome_dict = SAMTOOLS_DICT.out.dict
            .map { meta, dict -> [dict] }

        CHORD(
            ch_chord_input,
            chord_genome_fasta,
            chord_genome_fai,
            chord_genome_dict
        )
        ch_chord_prediction = CHORD.out.prediction
        ch_versions = ch_versions.mix(CHORD.out.versions)
    }

    //
    // MUTATIONALPATTERN: Mutational signature analysis
    //
    ch_mutational_profile = Channel.empty()

    if (!skip_mutationalpattern) {
        // Prepare genome parameter
        ch_mutationalpattern_genome = channel.of(
            [ [ id: 'hg38' ], 'BSgenome.Hsapiens.UCSC.hg38' ]
        )

        MUTATIONALPATTERN(
            ch_somatic_snv_vcf_gz,
            ch_mutationalpattern_genome,
            params.mutationalpattern_max_delta ?: 0.06
        )
        ch_mutational_signature = MUTATIONALPATTERN.out.mut_sig
        ch_versions = ch_versions.mix(MUTATIONALPATTERN.out.versions)
    }

    emit:
    chord_prediction      = ch_chord_prediction     // channel: [ meta, prediction.txt ]
    mutational_signature  = ch_mutational_profile   // channel: [ meta, signature ]
    versions              = ch_versions             // channel: [ versions.yml ]
}
