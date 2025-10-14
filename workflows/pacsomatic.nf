/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
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

include { PBTK_PBMERGE           } from '../modules/nf-core/pbtk/pbmerge/main'
include { PBMM2_ALIGN            } from '../modules/nf-core/pbmm2/align/main'
include { MOSDEPTH               } from '../modules/nf-core/mosdepth/main'
include	{ DEEPTOOLS_BAMCOVERAGE  } from	'../modules/nf-core/deeptools/bamcoverage/main'

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

    ch_processed = PBMM2_ALIGN.out.bam

    //
    // Pre-pare tumor-normal pairs for variant calling
    //

    // Split samples by status (0=normal, 1=tumor)
    ch_samples_by_patient = ch_processed
        .branch { meta, bam ->
            normal: meta.status == 0
                return [ meta.patient, meta, bam ]
            tumor: meta.status == 1
                return [ meta.patient, meta, bam ]
        }

    // Group normals and tumors by patient
    ch_normals = ch_samples_by_patient.normal
        .map { patient, meta, bam  ->
            [ patient, meta, bam ]
        }

    ch_tumors = ch_samples_by_patient.tumor
        .map { patient, meta, bam ->
            [ patient, meta, bam ]
        }

    // Generate tumor and normal pairs for each patient
    ch_tn_pairs = ch_tumors
        .combine(ch_normals, by: 0)
        .map { patient, tumor_meta, tumor_bam,
            normal_meta, normal_bam ->
            def pair_meta = [
                patient: patient,
                tumor_id: tumor_meta.sample,
                normal_id: normal_meta.sample,
                id: "${patient}_${tumor_meta.sample}_vs_${normal_meta.sample}"
            ]
            [ pair_meta, tumor_bam, normal_bam ]
        }
    // ch_tn_pairs.view()

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
