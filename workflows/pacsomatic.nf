/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNZIPFILES             } from '../modules/nf-core/unzipfiles/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pacsomatic_pipeline'
include { checkParameters        } from '../subworkflows/local/utils_pacsomatic_pipeline'
include { checkPathParameters    } from '../subworkflows/local/utils_pacsomatic_pipeline'

include { PREPARE_GENOME          } from '../subworkflows/local/prepare_genome'
include { BAM_SORT_STATS_SAMTOOLS } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

include { PBTK_PBMERGE           } from '../modules/local/pbtk/pbmerge/main'
include { PBMM2_ALIGN            } from '../modules/nf-core/pbmm2/align/main'

include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include	{ DEEPTOOLS_BAMCOVERAGE  } from	'../modules/nf-core/deeptools/bamcoverage/main'

include { DEEPSOMATIC            } from '../modules/nf-core/deepsomatic/main'
//include { MUTATIONALPATTERN    } from '../modules/local/mutationalpattern/main'
//include { ENSEMBLVEP_DOWNLOAD    } from '../modules/nf-core/ensemblvep/download/main'
//include { ENSEMBLVEP_VEP         } from '../modules/nf-core/ensemblvep/vep/main'

include { SEVERUS                } from '../modules/nf-core/severus/main'
//include { ANNOTSV_ANNOTSV        } from '../modules/nf-core/annotsv/annotsv/main'

include { CHORD                  } from '../modules/local/chord/main'

include { CNVKIT_BATCH           } from '../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL            } from '../modules/nf-core/cnvkit/call/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PACSOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // check parameters
    checkParameters()
    checkPathParameters()

    // init version and multiqc channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // Group BAMs by patient-sample and merge if multiple files exists
    //
    ch_grouped_bams = ch_samplesheet
        .map { meta, bam, pbi ->
            def patient_sample = "${meta.patient}_${meta.sample}"
            // Add id field for process naming
            def meta_with_id = meta + [id: patient_sample]
            [ patient_sample, meta_with_id, bam, pbi ]
        }
        .groupTuple(by: 0)
        .branch { _patient_sample, meta_list, bam_list, pbi_list ->
            // Take first meta since patient/sample/status should be identical
            def meta = meta_list[0]
            single: bam_list.size() == 1
                return [ meta, bam_list[0], pbi_list[0] ]
            multiple: bam_list.size() > 1
                return [ meta, bam_list, pbi_list ]
        }

    // Merge multiple BAMs per patient-sample
    PBTK_PBMERGE(
        ch_grouped_bams.multiple.map { meta, bams, _pbis ->
            [ meta, bams ]
        }
    )
    ch_versions = ch_versions.mix(PBTK_PBMERGE.out.versions.first())

    // Combine single BAMs and merged BAMs
    ch_input_sample = ch_grouped_bams.single
        .mix(
            PBTK_PBMERGE.out.bam.join(PBTK_PBMERGE.out.pbi, by: 0)
        )

    //
    // SUBWORKFLOW: Uncompress and prepare reference genomes files used by the pipeline
    //
    PREPARE_GENOME( params.fasta )
    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)

    ch_genome_fasta = PREPARE_GENOME.out.prepped_genome_fasta
    ch_genome_fai   = PREPARE_GENOME.out.genome_fai

    // Alignment with PacBio PBMM2
    ch_pbmm2_input = ch_input_sample
        .map { meta, bam, _pbi ->
            [ meta, bam ]
        }
    PBMM2_ALIGN ( ch_pbmm2_input, ch_genome_fasta )
    ch_versions = ch_versions.mix(PBMM2_ALIGN.out.versions.first())

    //
    /* YOUR ANALYSIS BEFORE PAIRING STARTS HERE */
    //
    
    BAM_SORT_STATS_SAMTOOLS(PBMM2_ALIGN.out.bam, ch_genome_fasta)

    // join the bam and index based on meta.id
    ch_ordered_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
    ch_ordered_bai = BAM_SORT_STATS_SAMTOOLS.out.bai

    ch_bam_bai = ch_ordered_bam.join(ch_ordered_bai, by: [0])
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // these files go for multiqc
    ch_ordered_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
    ch_ordered_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    ch_ordered_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats
    
    ch_mosdepth_multiqc_files = Channel.empty()
    if ( !params.skip_qc && !params.skip_mosdepth ) {
        ch_mosdepth_input = ch_bam_bai.map{ meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }
        MOSDEPTH (ch_mosdepth_input, ch_genome_fasta)

        ch_mosdepth_multiqc_files = ch_mosdepth_multiqc_files.mix(MOSDEPTH.out.global_txt)
        ch_mosdepth_multiqc_files = ch_mosdepth_multiqc_files.mix(MOSDEPTH.out.summary_txt)
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    }

    // generates a coverage track using deeptools/bamcoverage
    if ( !params.skip_qc && !params.skip_bamcoverage )  {
        DEEPTOOLS_BAMCOVERAGE(
            ch_bam_bai,
            ch_genome_fasta.map{it[1]},
            ch_genome_fai.map{it[1]},
            [[:], []])
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    }
    
    ch_processed = ch_bam_bai // PBMM2_ALIGN.out.bam

    //
    // Pre-pare tumor-normal pairs for variant calling
    //

    // Split samples by status (0=normal, 1=tumor)
    ch_samples_by_patient = ch_processed
        .branch { meta, bam, bai ->
            normal: meta.status == 0
                return [ meta.patient, meta, bam, bai ]
            tumor: meta.status == 1
                return [ meta.patient, meta, bam, bai ]
        }

    // Group normals and tumors by patient
    ch_normals = ch_samples_by_patient.normal
        .map { patient, meta, bam, bai  ->
            [ patient, meta, bam, bai ]
        }

    ch_tumors = ch_samples_by_patient.tumor
        .map { patient, meta, bam, bai ->
            [ patient, meta, bam, bai ]
        }

    // Generate tumor and normal pairs for each patient
    ch_tn_bam_pairs = ch_tumors
        .combine(ch_normals, by: 0)
        .map { patient, tumor_meta, tumor_bam, tumor_bam_bai,
            normal_meta, normal_bam, normal_bam_bai ->
            def pair_meta = [
                patient: patient,
                tumor_id: tumor_meta.sample,
                normal_id: normal_meta.sample,
                id: "${patient}_${tumor_meta.sample}_vs_${normal_meta.sample}"
            ]
            [ pair_meta, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ]
        }
    // ch_tn_bam_pairs.view()

    ch_somatic_snv_vcf_gz= Channel.empty()
    if (!params.skip_deepsomatic)  {
        ch_deepsomatic_tn_bam_pairs =  ch_tn_bam_pairs // as order  [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        ch_deepsomatic_interval     =  channel.of( [[:], []] )  // channel.of( [[:], "$projectDir/assets/dummy_file.bed"] )
        ch_deepsomatic_gzi          =  channel.of( [[:], []] )  // channel.of( [[:], "$projectDir/assets/dummy_file.gz"] )

        DEEPSOMATIC(ch_deepsomatic_tn_bam_pairs, ch_deepsomatic_interval, ch_genome_fasta, ch_genome_fai, ch_deepsomatic_gzi)
        ch_somatic_snv_vcf_gz=DEEPSOMATIC.out.vcf
        ch_versions          =ch_versions.mix(DEEPSOMATIC.out.versions)
    }

    /*
    if (!params.skip_deepsomatic && !params.skip_vep) {
        ch_vep_download=Channel.of( [[:], "${params.vep_assembly}", "${params.vep_species}", "${params.vep_cache_version}"]) // meta, assembly, species, cache_version
        ENSEMBLVEP_DOWNLOAD(ch_vep_download)
        vep_cache_path  =ENSEMBLVEP_DOWNLOAD.out.cache
                         .map { meta, path_prefix -> [path_prefix]
                         }

        ch_versions        = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

        ch_vep_somatic_snv_vcf_gz =ch_somatic_snv_vcf_gz.combine([[]]);

        ENSEMBLVEP_VEP(ch_vep_somatic_snv_vcf_gz, "${params.vep_assembly}", "${params.vep_species}", "${params.vep_cache_version}", vep_cache_path, ch_genome, [])
        ch_versions        = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    }
    */

    ch_somatic_sv_vcf = Channel.empty()
    if (!params.skip_severus) {
        ch_severus_phasing_vcf   = channel.of([[]])  // $projectDir/assets/dummy_file.vcf
        ch_severus_tn_bam_pairs  = ch_tn_bam_pairs
        .map {
               pair_meta,  normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ->
               [pair_meta, tumor_bam, tumor_bam_bai, normal_bam, normal_bam_bai]
            }
	.combine(ch_severus_phasing_vcf)  // as order of [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        //ch_severus_tn_bam_pairs.view()

        ch_severus_trf_bed	 = channel.of([[:],[]]) // channel.of([[:],"${params.Severus_trf_bed}"])   need a s3 bucket to store configurable bed
        SEVERUS(ch_severus_tn_bam_pairs, ch_severus_trf_bed)
        ch_somatic_sv_vcf = SEVERUS.out.somatic_vcf
                             .map {pair_meta, sv_vcf ->
                             [ pair_meta.id, pair_meta, sv_vcf ]
                             } // [pair_id, pair_meta, somatic_vcf]
        // ch_somatic_sv_vcf.view()
        ch_versions	  = ch_versions.mix(SEVERUS.out.versions)
    }
    
    if( !params.skip_deepsomatic && !params.skip_severus && !params.skip_chord) {
        UNZIPFILES (ch_somatic_snv_vcf_gz)
        ch_somatic_snv_vcf = UNZIPFILES.out.files
                             .map {pair_meta, snv_vcf ->
                             [pair_meta.id, pair_meta, snv_vcf]
                              }

        //ch_somatic_snv_vcf.view()
        //ch_somatic_sv_vcf.view()

        ch_chord_snv_sv_vcfs= ch_somatic_snv_vcf.combine(ch_somatic_sv_vcf, by :0)
                              .map { pair_id, pair_meta, snv_vcf, pair_meta1, sv_vcf ->
                               def chord_meta = [
                               patient: pair_meta.patient,
                               tumor_id: pair_meta.tumor_id,
                               normal_id: pair_meta.normal_id,
                               id: pair_meta.id,
                               sample_id: pair_meta.id 
                              ]
                              [chord_meta, snv_vcf, sv_vcf]
                              }

        ch_chord_snv_sv_vcfs.view()

        chord_genome_fasta  = ch_genome_fasta.map {meta, genome_fasta -> [genome_fasta] }
        chord_genome_fai    = ch_genome_fai.map { meta, genome_fai -> [genome_fai] }

        chord_genome_fasta.view()
        chord_genome_fai.view()

        CHORD(ch_chord_snv_sv_vcfs, chord_genome_fasta, chord_genome_fai, [])

        ch_versions        = ch_versions.mix(CHORD.out.versions)
    }
    
    /*
    if (!params.skip_mutationalpattern) {
       ch_mutationalpattern_vcf   = ch_somatic_snv_vcf_gz
       ch_mutationalpattern_genome= ch_genome_fasta
       MUTATIONALPATTERN(ch_mutationalpattern_vcf, ch_mutationalpattern_genome, "${params.mutationalpattern_max_delta}")
       ch_versions        = ch_versions.mix(MUTATIONALPATTERN.out.versions)
    }
    */

    if (!params.skip_cnvkit) {
       //  First step as  CNVKit Batch
       ch_cnvkit_tn_bam_pairs = ch_tn_bam_pairs
       .map {
           pair_meta,  normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ->
           [pair_meta, tumor_bam, normal_bam]
       }

       // ch_cnvkit_tn_bam_pairs.view()

       ch_cnvkit_targets  = channel.of([[:],[]]) // need a s3 bucket to store configurable bed
       ch_cnvkit_reference= channel.of([[:],[]])
       CNVKIT_BATCH(ch_cnvkit_tn_bam_pairs, ch_genome_fasta, ch_genome_fai, ch_cnvkit_targets, ch_cnvkit_reference, [:] )
       ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)

       //should add optional refinement steps here, e.g Merge two clari3 called vcfs of a pair into a whole vcf. 
       ch_merged_germline_vcf= channel.of([[]])

       // Last step as CNVKit Call

       ch_cnvkit_cns = CNVKIT_BATCH.out.cns   //[meta, [binter_cns, call_cns, cns] ]
       .combine(ch_merged_germline_vcf)       //[meta, [binter_cns, call_cns, cns], merge_vcf]
       .map {cnvkit_meta, cns, merge_vcf ->
       [cnvkit_meta, cns[1], merge_vcf]       //[meta, call_cns, merge_vcf]
       }

       // ch_cnvkit_cns.view()

       CNVKIT_CALL(ch_cnvkit_cns)
       ch_versions = ch_versions.mix(CNVKIT_CALL.out.versions)

    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'pacsomatic_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    if ( !params.skip_qc && !params.skip_multiqc ) {

        //
        // MODULE: MultiQC
        //
        ch_multiqc_config        = Channel.fromPath(
            "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ?
            Channel.fromPath(params.multiqc_config, checkIfExists: true) :
            Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo ?
            Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
            Channel.empty()

        summary_params      = paramsSummaryMap(
            workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
            file(params.multiqc_methods_description, checkIfExists: true) :
            file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
        ch_methods_description                = Channel.value(
            methodsDescriptionText(ch_multiqc_custom_methods_description))

        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        ch_multiqc_files = ch_multiqc_files.mix(
            ch_methods_description.collectFile(
                name: 'methods_description_mqc.yaml',
                sort: true
            )
        )

        // post-alignment qc files

        ch_multiqc_files = ch_multiqc_files.mix(ch_ordered_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_ordered_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_ordered_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_multiqc_files.collect{it[1]}.ifEmpty([]))

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            [],
            []
        )

        ch_multiqc_report = MULTIQC.out.report
        ch_versions = ch_versions.mix(MULTIQC.out.versions)

    }

    emit:
    multiqc_report = ch_multiqc_report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
