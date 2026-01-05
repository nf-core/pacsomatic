/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText     } from '../subworkflows/local/utils_nfcore_pacsomatic_pipeline'
include { checkParameters            } from '../subworkflows/local/utils_pacsomatic_pipeline'
include { checkPathParameters        } from '../subworkflows/local/utils_pacsomatic_pipeline'

// Reference genome preparation
include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'

// Alignment and QC
include { PBTK_PBMERGE               } from '../modules/nf-core/pbtk/pbmerge/main'
include { PBMM2_ALIGN                } from '../modules/nf-core/pbmm2/align/main'
include { BAM_SORT_STATS_SAMTOOLS    } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { MOSDEPTH                   } from '../modules/nf-core/mosdepth/main'
include { DEEPTOOLS_BAMCOVERAGE      } from '../modules/nf-core/deeptools/bamcoverage/main'

// Analysis subworkflows
include { GERMLINE_CALLING_PHASING   } from '../subworkflows/local/germline_calling_phasing/main'
include { SOMATIC_SNV_INDEL          } from '../subworkflows/local/somatic_snv_indel/main'
include { SOMATIC_SV                 } from '../subworkflows/local/somatic_sv/main'
include { CNV_CALLING                } from '../subworkflows/local/cnv_calling/main'
include { TUMOR_CLONALITY            } from '../subworkflows/local/tumor_clonality/main'
include { METHYLATION_ANALYSIS       } from '../subworkflows/local/methylation_analysis/main'
include { SIGNATURE_ANALYSIS         } from '../subworkflows/local/signature_analysis/main'

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
    ch_multiqc_report = Channel.empty()

    //
    // initialize file channels for subworkflows
    //
    // SV-related channels
    ch_severus_trf_bed = params.severus_trf_bed
        ? Channel.value([[:], file(params.severus_trf_bed, checkIfExists: true)])
        : Channel.value([[:], []])

    ch_svpack_ctrl_vcf = params.svpack_ctrl_vcf
        ? Channel.value([[:], file(params.svpack_ctrl_vcf, checkIfExists: true)])
        : Channel.value([[:], []])

    ch_svpack_ref_gff = params.svpack_ref_gff
        ? Channel.value([[:], file(params.svpack_ref_gff, checkIfExists: true)])
        : Channel.value([[:], []])

    ch_annotsv_cache = params.annotsv_cache
        ? channel.value([[:], file(params.annotsv_cache, checkIfExists: true)])
        : Channel.empty()

    // CNV channels
    ch_cnv_target_bed = params.cnv_target_bed
        ? Channel.value([[:], file(params.cnv_target_bed, checkIfExists: true)])
        : Channel.value([[:], []])

    ch_cnv_reference = params.cnv_reference
        ? Channel.value([[:], file(params.cnv_reference, checkIfExists: true)])
        : Channel.value([[:], []])

    ch_cnv_germline_vcf = params.cnv_germline_vcf
        ? Channel.value(file(params.cnv_germline_vcf, checkIfExists: true))
        : Channel.value([[]])

    // Tumor clonality channels
    ch_target_regions_bed = params.target_regions_bed
        ? Channel.value(file(params.target_regions_bed, checkIfExists: true))
        : Channel.value([])

    ch_heterozygous_sites = params.heterozygous_sites
        ? Channel.value(file(params.heterozygous_sites, checkIfExists: true))
        : Channel.value([])

    ch_diploid_regions = params.diploid_regions
        ? Channel.value(file(params.diploid_regions, checkIfExists: true))
        : Channel.value([])

    ch_target_region_normalisation = params.target_region_normalisation
        ? Channel.value(file(params.target_region_normalisation, checkIfExists: true))
        : Channel.value([])

    ch_gc_profile = params.gc_profile
        ? Channel.value(file(params.gc_profile, checkIfExists: true))
        : Channel.value([])

    // Variant calling channels
    ch_known_hotspots_somatic = params.known_hotspots_somatic
        ? Channel.value(file(params.known_hotspots_somatic, checkIfExists: true))
        : Channel.value([])

    ch_known_hotspots_germline = params.known_hotspots_germline
        ? Channel.value(file(params.known_hotspots_germline, checkIfExists: true))
        : Channel.value([])

    ch_driver_gene_panel = params.driver_gene_panel
        ? Channel.value(file(params.driver_gene_panel, checkIfExists: true))
        : Channel.value([])

    ch_ensembl_data_dir = params.ensembl_data_dir
        ? Channel.value(file(params.ensembl_data_dir, checkIfExists: true, type: 'dir'))
        : Channel.value([])

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
    // SUBWORKFLOW: Post-alignment processing (sorting, indexing, stats)
    //
    BAM_SORT_STATS_SAMTOOLS(PBMM2_ALIGN.out.bam, ch_genome_fasta)
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // Join BAM and BAI
    ch_bam_bai = BAM_SORT_STATS_SAMTOOLS.out.bam
        .join(BAM_SORT_STATS_SAMTOOLS.out.bai, by: [0])

    // Collect QC files for MultiQC
    ch_ordered_stats    = BAM_SORT_STATS_SAMTOOLS.out.stats
    ch_ordered_flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat
    ch_ordered_idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats

    //
    // MODULE: Coverage analysis with MOSDEPTH (optional)
    //
    ch_mosdepth_multiqc_files = Channel.empty()
    if (!params.skip_qc && !params.skip_mosdepth) {
        ch_mosdepth_input = ch_bam_bai.map { meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }
        MOSDEPTH(ch_mosdepth_input, ch_genome_fasta)

        ch_mosdepth_multiqc_files = ch_mosdepth_multiqc_files.mix(MOSDEPTH.out.global_txt)
        ch_mosdepth_multiqc_files = ch_mosdepth_multiqc_files.mix(MOSDEPTH.out.summary_txt)
        ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first())
    }

    //
    // MODULE: Generate coverage tracks with DEEPTOOLS (optional)
    //
    if (!params.skip_qc && !params.skip_bamcoverage) {
        DEEPTOOLS_BAMCOVERAGE(
            ch_bam_bai,
            ch_genome_fasta.map { it[1] },
            ch_genome_fai.map { it[1] },
            [[:], []]
        )
        ch_versions = ch_versions.mix(DEEPTOOLS_BAMCOVERAGE.out.versions.first())
    }

    //
    // SUBWORKFLOW: Germline variant calling and phasing
    //
    GERMLINE_CALLING_PHASING(
        ch_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        params.skip_hiphase
    )
    ch_versions = ch_versions.mix(GERMLINE_CALLING_PHASING.out.versions)

    ch_vcf_tbi = GERMLINE_CALLING_PHASING.out.vcf
    ch_phased_normal_bam = GERMLINE_CALLING_PHASING.out.phased_bam

    //
    // Create tumor-normal pairs for somatic analysis
    // Split samples by status and create paired channels
    //
    ch_samples_by_patient = ch_bam_bai
        .branch { meta, bam, bai ->
            normal: meta.status == 0
                return [ meta.patient, meta, bam, bai ]
            tumor: meta.status == 1
                return [ meta.patient, meta, bam, bai ]
        }

    // Prepare tumor and normal BAM channels
    ch_bam_normals_tmp = ch_samples_by_patient.normal
    ch_bam_tumors_tmp  = ch_samples_by_patient.tumor

    ch_bam_tumors = ch_bam_tumors_tmp
        .map { _patient, meta, bam, bai ->
            [ meta, bam, bai ]
        }

    // Create tumor-normal paired BAM channel
    // Structure: [ pair_meta, normal_bam, normal_bai, tumor_bam, tumor_bai ]
    ch_tn_bam_pairs = ch_bam_tumors_tmp
        .combine(ch_bam_normals_tmp, by: 0)
        .map { patient, tumor_meta, tumor_bam, tumor_bam_bai,
            normal_meta, normal_bam, normal_bam_bai ->
            def pair_meta = [
                patient:   patient,
                tumor_id:  tumor_meta.sample,
                normal_id: normal_meta.sample,
                id:        "${patient}_${tumor_meta.sample}_vs_${normal_meta.sample}"
            ]
            [ pair_meta, normal_bam, normal_bam_bai, tumor_bam, tumor_bam_bai ]
        }

    // Prepare tumor and normal VCF channels (germline variants)
    ch_vcf_by_patient = ch_vcf_tbi
        .branch { meta, vcf, tbi ->
            normal: meta.status == 0
                return [ meta.patient, meta, vcf, tbi ]
            tumor: meta.status == 1
                return [ meta.patient, meta, vcf, tbi ]
        }

    ch_vcf_normals_tmp = ch_vcf_by_patient.normal
    ch_vcf_tumors_tmp  = ch_vcf_by_patient.tumor

    ch_vcf_normals = ch_vcf_normals_tmp
        .map { patient, meta, vcf, tbi ->
            [ meta, vcf, tbi ]
        }

    ch_vcf_tumors = ch_vcf_tumors_tmp
        .map { patient, meta, vcf, tbi ->
            [ meta, vcf, tbi ]
        }

    //
    // SUBWORKFLOW: Somatic SNV/indel calling, annotation, and phasing
    //
    ch_somatic_snv_vcf_gz = Channel.empty()
    ch_somatic_phased_bam = Channel.empty()

    if (!params.skip_deepsomatic) {
        SOMATIC_SNV_INDEL(
            ch_tn_bam_pairs,
            ch_genome_fasta,
            ch_genome_fai,
            ch_vcf_tumors,
            ch_bam_tumors,
            params.skip_ensemblvep,
            params.skip_vep_download,
            params.skip_somatic_hiphase
        )
        ch_versions = ch_versions.mix(SOMATIC_SNV_INDEL.out.versions)

        ch_somatic_snv_vcf_gz  = SOMATIC_SNV_INDEL.out.vcf
        ch_somatic_phased_bam  = SOMATIC_SNV_INDEL.out.phased_bam
    }

    //
    // SUBWORKFLOW: Somatic structural variant calling and annotation
    //
    ch_somatic_sv_vcf = Channel.empty()

    if (!params.skip_sv) {
        SOMATIC_SV(
            ch_tn_bam_pairs,
            ch_severus_trf_bed,
            ch_svpack_ctrl_vcf,
            ch_svpack_ref_gff,
            ch_annotsv_cache,
            params.skip_svpack,
            params.skip_annotsv,
            params.skip_annotsv_install
        )

        ch_somatic_sv_vcf = SOMATIC_SV.out.vcf
        ch_versions = ch_versions.mix(SOMATIC_SV.out.versions)
    }

    //
    // SUBWORKFLOW: Copy number variant calling
    //
    if (!params.skip_cnvkit) {
        CNV_CALLING(
            ch_tn_bam_pairs,
            ch_genome_fasta,
            ch_genome_fai,
            ch_cnv_target_bed,
            ch_cnv_reference,
            ch_cnv_germline_vcf,
        )
        ch_versions = ch_versions.mix(CNV_CALLING.out.versions)
    }

    //
    // SUBWORKFLOW: Tumor clonality analysis (AMBER, COBALT, PURPLE)
    //
    if (!params.skip_tumor_clonality) {
        TUMOR_CLONALITY(
            ch_tn_bam_pairs,
            ch_genome_fasta,
            ch_genome_fai,
            ch_heterozygous_sites,
            ch_target_regions_bed,
            ch_gc_profile,
            ch_diploid_regions,
            ch_target_region_normalisation,
            ch_known_hotspots_somatic,
            ch_known_hotspots_germline,
            ch_driver_gene_panel,
            ch_ensembl_data_dir
        )
        ch_versions = ch_versions.mix(TUMOR_CLONALITY.out.versions)
    }

    //
    // SUBWORKFLOW: Methylation analysis and DMR detection
    //
    // Methylation analysis requires phased BAMs from both germline and somatic phasing
    if (!params.skip_pbcpgtools && !params.skip_hiphase && !params.skip_somatic_hiphase) {
        METHYLATION_ANALYSIS(
            ch_phased_normal_bam,
            ch_somatic_phased_bam,
            params.skip_dmr,
            params.skip_dmr_anno
        )
        ch_versions = ch_versions.mix(METHYLATION_ANALYSIS.out.versions)
    }

    //
    // SUBWORKFLOW: Mutational signature analysis (CHORD and MUTATIONALPATTERN)
    //
    if ((!params.skip_chord || !params.skip_mutationalpattern) &&
        !params.skip_deepsomatic && !params.skip_sv) {
        SIGNATURE_ANALYSIS(
            ch_somatic_snv_vcf_gz,
            ch_somatic_sv_vcf,
            ch_genome_fasta,
            ch_genome_fai,
            params.skip_chord,
            params.skip_mutationalpattern
        )
        ch_versions = ch_versions.mix(SIGNATURE_ANALYSIS.out.versions)
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
