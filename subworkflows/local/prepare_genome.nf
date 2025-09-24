//
// Uncompress and prepare reference genome files
//

include { PIGZ_UNCOMPRESS as GUNZIP_GENOME_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { UNZIPFILES as UNZIP_GENOME_FASTA       } from '../../modules/nf-core/unzipfiles/main'


workflow PREPARE_GENOME {

    take:
    genome_fasta  // str: path/to/genome.fasta

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: Unzip Genome Fasta
    //
    ch_genome_fasta = Channel.empty()

    if (genome_fasta) {
        if (genome_fasta.endsWith('.gz')) {
            GUNZIP_GENOME_FASTA ([ [:], file(genome_fasta, checkIfExists: true) ])

            ch_genome_fasta = GUNZIP_GENOME_FASTA.out.file
            ch_versions = ch_versions.mix(GUNZIP_GENOME_FASTA.out.versions)

        } else if (genome_fasta.endsWith('.zip')) {
            UNZIP_GENOME_FASTA ([ [:], file(genome_fasta, checkIfExists: true) ])
        } else {
            ch_genome_fasta = [ [:], file(genome_fasta, checkIfExists: true) ]
        }
    }

    emit:
    prepped_genome_fasta = ch_genome_fasta
    versions             = ch_versions
}
