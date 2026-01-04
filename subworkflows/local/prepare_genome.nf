//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_GENOME_FASTA    } from '../../modules/nf-core/gunzip/main'
include { UNZIPFILES as UNZIP_GENOME_FASTA } from '../../modules/nf-core/unzipfiles/main'
include { SAMTOOLS_FAIDX as GENOME_FAIDX   } from '../../modules/nf-core/samtools/faidx/main'


workflow PREPARE_GENOME {

    take:
    genome_fasta  // str: path/to/genome.fasta

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Unzip Genome Fasta
    //
    ch_genome_fasta = Channel.empty()
    ch_genome_fai = Channel.empty()

    if (genome_fasta) {
        if (genome_fasta.endsWith('.gz')) {
            GUNZIP_GENOME_FASTA ([ [:], file(genome_fasta, checkIfExists: true) ])

            ch_genome_fasta = GUNZIP_GENOME_FASTA.out.file
            ch_versions = ch_versions.mix(GUNZIP_GENOME_FASTA.out.versions)

        } else if (genome_fasta.endsWith('.zip')) {
            UNZIP_GENOME_FASTA ([ [:], file(genome_fasta, checkIfExists: true) ])

            ch_genome_fasta = UNZIP_GENOME_FASTA.out.files
            ch_versions = ch_versions.mix(UNZIP_GENOME_FASTA.out.versions)

        } else {
            ch_genome_fasta = [ [:], file(genome_fasta, checkIfExists: true) ]
        }

        //
        // MODULE: Index the genome fasta
        //
        GENOME_FAIDX( ch_genome_fasta, [ [:], "$projectDir/assets/dummy_file.txt" ], [] )
        ch_genome_fai = GENOME_FAIDX.out.fai
        ch_versions = ch_versions.mix(GENOME_FAIDX.out.versions)
    }

    emit:
    prepped_genome_fasta = ch_genome_fasta
    genome_fai           = ch_genome_fai
    versions             = ch_versions
}
