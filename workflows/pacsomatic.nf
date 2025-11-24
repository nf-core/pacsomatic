/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNZIPFILES                 } from '../modules/nf-core/unzipfiles/main'
include { FASTQC                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                    } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc       } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText     } from '../subworkflows/local/utils_nfcore_pacsomatic_pipeline'
include { checkParameters            } from '../subworkflows/local/utils_pacsomatic_pipeline'
include { checkPathParameters        } from '../subworkflows/local/utils_pacsomatic_pipeline'

include { PREPARE_GENOME             } from '../subworkflows/local/prepare_genome'
include { BAM_SORT_STATS_SAMTOOLS    } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

include { PBTK_PBMERGE               } from '../modules/nf-core/pbtk/pbmerge/main'
include { PBMM2_ALIGN                } from '../modules/nf-core/pbmm2/align/main'
include { MOSDEPTH                   } from '../modules/nf-core/mosdepth/main'
include	{ DEEPTOOLS_BAMCOVERAGE      } from '../modules/nf-core/deeptools/bamcoverage/main'
include { CLAIR3                     } from '../modules/nf-core/clair3/main'

include { DEEPSOMATIC                } from '../modules/nf-core/deepsomatic/main'
include { MUTATIONALPATTERN          } from '../modules/local/mutationalpattern/main'
include { ENSEMBLVEP_DOWNLOAD        } from '../modules/nf-core/ensemblvep/download/main'
include { ENSEMBLVEP_VEP             } from '../modules/nf-core/ensemblvep/vep/main'
include { AMBER                      } from '../modules/local/amber/main'
include { COBALT                     } from '../modules/local/cobalt/run/main'
//include { COBALT_PANEL_NORMALISATION } from '../modules/local/cobalt/panel_normalisation/main'
include { PURPLE                     } from '../modules/local/purple/main' 

include { SEVERUS                    } from '../modules/nf-core/severus/main'
include { SVPACK_ANNOTATE            } from '../subworkflows/local/svpack_annotate/main'
include { ANNOTSV_INSTALLANNOTATIONS } from '../modules/nf-core/annotsv/installannotations/main'
include { ANNOTSV_ANNOTSV            } from '../modules/nf-core/annotsv/annotsv/main'
include { TABIX_BGZIPTABIX as TABIX_SV_VCF } from '../modules/nf-core/tabix/bgziptabix/main'

include { SAMTOOLS_DICT              } from '../modules/nf-core/samtools/dict/main'
include { CHORD                      } from '../modules/local/chord/main'

include { CNVKIT_BATCH               } from '../modules/nf-core/cnvkit/batch/main'
include { CNVKIT_CALL                } from '../modules/nf-core/cnvkit/call/main'

include { HIPHASE                    } from '../modules/nf-core/hiphase/main'
include { HIPHASE_SOMATIC            } from '../modules/local/hiphase_somatic/main'
// include	{ HIPHASE as HIPHASE_SOMATIC } from '../modules/nf-core/hiphase/main'

include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES as PBCPGTOOLS_NORMAL} from '../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
include { PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES as PBCPGTOOLS_TUMOR } from '../modules/nf-core/pbcpgtools/alignedbamtocpgscores/main'
 
include { DSS_DMR                    } from '../modules/local/dss_dmr/main'
include { ANNOTATR_DMR               } from '../modules/local/annotatr_dmr/main'

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
    // Germline variant calling with CLAIR3
    //
    if (params.clair3_model && params.clair3_model_path) {
        log.error "Two models specified ${params.clair3_model} and ${params.clair3_model_path}, specify one of them."
        exit 1
    }
    if (!params.clair3_model && !params.clair3_model_path) {
        log.error "No clair3 model is specified, use option params.clair3_model or params.clair3_model_path."
        exit 1
    }
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
    
    ch_vcf_tbi=CLAIR3.out.vcf.join(CLAIR3.out.tbi)    
   
    //
    // Pre-pare tumor-normal pairs for variant calling.
    // Generating the tumor , normal and tumot-normal bam channels
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
    ch_bam_normals_tmp = ch_samples_by_patient.normal

    ch_bam_normals     = ch_bam_normals_tmp
        .map { patient, meta, bam, bai  ->
            [ meta, bam, bai ]
        }

    
    ch_bam_tumors_tmp = ch_samples_by_patient.tumor

    ch_bam_tumors     = ch_bam_tumors_tmp
        .map { patient, meta, bam, bai ->
            [ meta, bam, bai ]
        }

    // Generate tumor and normal pairs for each patient
    ch_tn_bam_pairs = ch_bam_tumors_tmp
        .combine(ch_bam_normals_tmp, by: 0)
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
    
    // Generating the tumor, normal, and tumor-normal vcf channels
    ch_vcf_by_pateint = ch_vcf_tbi
                       .branch { meta, vcf, tbi ->
                           normal: meta.status == 0
                               return [ meta.patient, meta, vcf, tbi ]
                           tumor: meta.status == 1
                               return [ meta.patient, meta, vcf, tbi ]
                        }

    ch_vcf_normals_tmp = ch_vcf_by_pateint.normal

    ch_vcf_normals = ch_vcf_normals_tmp 
                    .map { patient, meta, vcf, tbi ->
                           [ meta, vcf, tbi ]
                    }

    ch_vcf_tumors_tmp = ch_vcf_by_pateint.tumor

    ch_vcf_tumors     = ch_vcf_tumors_tmp 
                    .map { patient, meta, vcf, tbi ->
                           [ meta, vcf, tbi ]
                    }
 
    ch_tn_vcf_pair = ch_vcf_tumors_tmp
                     .combine(ch_vcf_normals_tmp, by: 0)
                     .map{ patient, tumor_meta, tumor_vcf, tumor_vcf_tbi,
                           normal_meta, normal_vcf, normal_vcf_tbi ->
                       def pair_meta = [
                       patient: patient,
                       tumor_id: tumor_meta.sample,
                       normal_id: normal_meta.sample,
                       id: "${patient}_${tumor_meta.sample}_vs_${normal_meta.sample}"
                       ]
                       [ pair_meta, normal_vcf, normal_vcf_tbi, tumor_vcf, tumor_vcf_tbi ]
                    }
    //
    //  Germline HiPhase for normal channel
    //
    ch_hiphase_out_bam_bai= Channel.empty()
    if (!params.skip_hiphase) {
        ch_hiphase_vcf     =  ch_vcf_normals     //  ch_vcf_tbi //= CLAIR3.out.vcf.join(CLAIR3.out.tbi)
                             .map { meta, vcf, tbi ->
                             [meta.id, meta, vcf, tbi]
                             }

        ch_hiphase_bam_bai = ch_bam_normals     // ch_bam_bai
                            .map { meta, bam, bai ->
                             [meta.id, meta, bam, bai]
                             }        
        
        ch_hiphase_combine = ch_hiphase_vcf.combine(ch_hiphase_bam_bai, by: [0])
                             .multiMap { meta_id, meta, vcf, tbi, meta2, bam, bai ->     
                               vcf_tbi: [meta, vcf, tbi]
                               bam_bai: [meta2, bam, bai]
                              }
         
        HIPHASE (ch_hiphase_combine.vcf_tbi, ch_hiphase_combine.bam_bai, ch_genome_fasta)
        ch_hiphase_out_bam_bai= HIPHASE.out.bam.join(HIPHASE.out.bai)
        ch_versions = ch_versions.mix(HIPHASE.out.versions.first())
        
    } 
    
    //
    //  TUMOR CLONALITY using hfmtools: amber, cobalt and purple   
    //
    if (!params.skip_tumor_clonality) {
        ch_amber =  ch_tn_bam_pairs  
                    .map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
                      [meta, tumor_bam, normal_bam, [], tumor_bai, normal_bai, [] ] // [ meta, tumor_bam, normal_bam, donor_bam, tumor_bai, normal_bai, donor_bai ]
                    }

       // ch_amber.view()

       // AMBER(ch_amber, 'V38', params.heterozygous_sites, params.target_regions_bed, params.tumor_min_depth)
       AMBER(ch_amber, 'V38', params.heterozygous_sites, params.target_regions_bed, [])
       //  AMBER(ch_amber, 'V38', params.heterozygous_sites, [], [])   //work
       ch_versions          =ch_versions.mix(AMBER.out.versions)
       ch_aber_dir         =AMBER.out.amber_dir
                            .map { pair_meta, amber_dir ->
                             [pair_meta.id, pair_meta, amber_dir] 
                            } 
      
       ch_cobalt = ch_tn_bam_pairs
                   .map { meta, normal_bam, normal_bai, tumor_bam, tumor_bai ->
                      [ meta, tumor_bam, normal_bam, tumor_bai, normal_bai ] 
                   }
       
       COBALT(ch_cobalt, params.gc_profile, params.diploid_regions, params.target_region_normalisation, [:])
       // COBALT(ch_cobalt, params.gc_profile, [], [],[:]) //work
       ch_versions          =ch_versions.mix(COBALT.out.versions)
       ch_cobalt_dir        =COBALT.out.cobalt_dir
                             .map { pair_meta, cobalt_dir ->
                               [pair_meta.id, pair_meta, cobalt_dir ]
                             }  
       
       //  sv_hard_vcf= []
       // sv_hard_vcf_index =[]
       // sv_soft_vcf=[]  
       // sv_soft_vcf_index=[]
       // smlv_tumor_vcf=[]
       // smlv_normal_vcf=[]

       ch_purple_amber_cobalt = AMBER.out.amber_dir.join(COBALT.out.cobalt_dir)
                                .map { meta, amber_dir, cobalt_dir ->
                                [ meta, amber_dir, cobalt_dir, [], [], [], [], [], [] ]
                                }

       PURPLE(ch_purple_amber_cobalt, ch_genome_fasta, ch_genome_fai, [[:],[]], '38', params.gc_profile, params.known_hotspots_somatic, params.known_hotspots_germline, params.driver_gene_panel, params.ensembl_data_dir, []) 
      
       // ch_cobalt_panel_normalisation=ch_amber_dir.combine(ch_cobalt_dir, by:[0])
       //                              .map { pair_meta_id, amber_pair_meta, amber_dir, cobalt_pair_meta, cobalt_dir ->
       //                               [amber_dir, cobalt_dir]
       //                               }
       
       //  COBALT_PANEL_NORMALISATION(ch_cobalt_panel_normalisation, 'V38', params.gc_profile, params.target_regions_bed)
       // ch_versions          =ch_versions.mix(COBALT_PANEL_NORMALISATION.out.versions) 

    }

    //
    // DEEPSOMATIC for Somatic SNV_INDEL calling
    //
    ch_somatic_snv_vcf_gz= Channel.empty()
    ch_somatic_snv_vcf_tbi= Channel.empty()
    if (!params.skip_deepsomatic)  {
        ch_deepsomatic_tn_bam_pairs =  ch_tn_bam_pairs // as order  [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
        ch_deepsomatic_interval     =  channel.of( [[:], []] )  // channel.of( [[:], "$projectDir/assets/dummy_file.bed"] )
        ch_deepsomatic_gzi          =  channel.of( [[:], []] )  // channel.of( [[:], "$projectDir/assets/dummy_file.gz"] )
        
        // ch_deepsomatic_tn_bam_pairs.view()

        DEEPSOMATIC(ch_deepsomatic_tn_bam_pairs, ch_deepsomatic_interval, ch_genome_fasta, ch_genome_fai, ch_deepsomatic_gzi)
        ch_somatic_snv_vcf_gz=DEEPSOMATIC.out.vcf
        ch_somatic_snv_vcf_tbi= DEEPSOMATIC.out.vcf_tbi
        ch_versions          =ch_versions.mix(DEEPSOMATIC.out.versions)
    }

    // 
    // VEP for somatic SNV annotation
    //
    if (!params.skip_deepsomatic && !params.skip_ensemblvep) {
        
        vep_cache_path="${params.vep_cache}"
        if ( !params.skip_vep_download ) {        
            ch_vep_download=Channel.of( [[:], "${params.vep_assembly}", "${params.vep_species}", "${params.vep_cache_version}"]) // meta, assembly, species, cache_version

            ENSEMBLVEP_DOWNLOAD(ch_vep_download)
            ch_versions        = ch_versions.mix(ENSEMBLVEP_DOWNLOAD.out.versions)

            vep_cache_path  =ENSEMBLVEP_DOWNLOAD.out.cache
                            .map { meta, path_prefix -> [path_prefix]
                             }
        }

        ch_vep_somatic_snv_vcf_gz =ch_somatic_snv_vcf_gz
                                   .map { meta, vcf_gz ->
                                   [ meta, vcf_gz, []]
                                   }

       // ch_vep_somatic_snv_vcf_gz.view()

        ENSEMBLVEP_VEP(ch_vep_somatic_snv_vcf_gz, "${params.vep_assembly}", "${params.vep_species}", "${params.vep_cache_version}", vep_cache_path, ch_genome_fasta, [])
        ch_versions        = ch_versions.mix(ENSEMBLVEP_VEP.out.versions)
    }

   //
   // SOMATIC Hiphasing for tumor channel   
   //
   ch_somatic_hiphase_out_bam_bai= Channel.empty()
   if ( !params.skip_deepsomatic && !params.skip_somatic_hiphase) {
       ch_tumor_hiphase_vcf     = ch_vcf_tumors     //  ch_vcf_tbi //= CLAIR3.out.vcf.join(CLAIR3.out.tbi)
                             .map { meta, vcf, tbi ->
                             [meta.id, meta, vcf, tbi]
                             }
      

       ch_somatic_hiphasing_vcf = ch_somatic_snv_vcf_gz.join(ch_somatic_snv_vcf_tbi)
                                  .map { pair_meta, snv_vcf, vcf_tbi ->
                                  def patient_tumor_id= "${pair_meta.patient}_${pair_meta.tumor_id}"
                                  [ patient_tumor_id, pair_meta, snv_vcf, vcf_tbi]
                                  }
        
       ch_somatic_hiphasing_bam_bai= ch_bam_tumors
                                     .map { meta, bam, bai ->
                                      [meta.id, meta, bam, bai]
                                    }
      
       ch_somatic_hiphasing_combine= ch_tumor_hiphase_vcf.combine(ch_somatic_hiphasing_bam_bai, by: 0 ).combine(ch_somatic_hiphasing_vcf, by: 0 )
                                     .multiMap {meta_id,  meta, vcf, tbi, meta2, bam, bai, meta3, somatic_vcf, somatic_tbi ->
                                     bam_bai: [meta2, bam, bai]
                                     vcf_tbi: [meta, vcf, tbi]
                                     somatic_vcf_tbi: [meta3, somatic_vcf, somatic_tbi]
                                    }

      HIPHASE_SOMATIC(ch_somatic_hiphasing_combine.vcf_tbi, ch_somatic_hiphasing_combine.bam_bai, ch_genome_fasta, ch_somatic_hiphasing_combine.somatic_vcf_tbi) 
      ch_somatic_hiphase_out_bam_bai = HIPHASE_SOMATIC.out.bam.join(HIPHASE_SOMATIC.out.bai)
      ch_versions                    = ch_versions.mix(HIPHASE_SOMATIC.out.versions)      
   }
   
   //
   //   pb_cpg methylation calling
   //
   if (!params.skip_pbcpgtools) {

       PBCPGTOOLS_NORMAL(ch_hiphase_out_bam_bai)
       ch_normal_cpg_bed = PBCPGTOOLS_NORMAL.out.combined_bed
                           .map { meta, cpg_bed ->
                           [meta.patient,  meta, cpg_bed]
                           }

       ch_versions       = ch_versions.mix(PBCPGTOOLS_NORMAL.out.versions.first())       
              
       
       PBCPGTOOLS_TUMOR(ch_somatic_hiphase_out_bam_bai)
       ch_tumor_cpg_bed = PBCPGTOOLS_TUMOR.out.combined_bed
                          .map{ meta,	cpg_bed	->
       	       	       	   [meta.patient, meta, cpg_bed]
                          }

       ch_versions      = ch_versions.mix(PBCPGTOOLS_TUMOR.out.versions.first())

       //
       //  Differential methylation region detection and annotation
       //
       if ( !params.skip_dmr && !params.skip_hiphase && !params.skip_deepsomatic && !params.skip_somatic_hiphase ) {
           
           // Differential methylation region detection
           ch_tn_pair_dmr_bed = ch_tumor_cpg_bed
           .combine(ch_normal_cpg_bed, by:[0])
           .map { patient, meta, tumor_bed, meta2, normal_bed ->
                [meta, tumor_bed, normal_bed]
           }

           DSS_DMR(ch_tn_pair_dmr_bed)
           ch_versions = ch_versions.mix(DSS_DMR.out.versions.first())
           
           // Annotation of detected differential methylation regions
           if ( !params.skip_dmr_anno )
           {           
              ch_dss_dmr_tsv= DSS_DMR.out.dmr
              ANNOTATR_DMR(ch_dss_dmr_tsv)
              ch_versions = ch_versions.mix(ANNOTATR_DMR.out.versions.first())
           }

       } 
   }
  
   //
   // SEVERUS for Somatic SV calling
   //
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
        // ch_somatic_sv_vcf.view()
        ch_versions	  = ch_versions.mix(SEVERUS.out.versions)
   }
   
   // 
   // SV Pack filtering 
   //
   ch_SV_annotsv_vcf =ch_somatic_sv_vcf
   if ( !params.skip_severus && !params.skip_svpack)
   { 
       ch_SVPack_vcf= ch_somatic_sv_vcf

       ch_SVPack_control_vcf = channel.of([[:],"${params.SVPack_control_vcf}"])
       ch_SVPack_ref_gff= channel.of([[:],"${params.SVPack_ref_gff}"])    

       SVPACK_ANNOTATE( ch_SVPack_vcf, ch_SVPack_control_vcf, ch_SVPack_ref_gff )
       ch_versions       = ch_versions.mix(SVPACK_ANNOTATE.out.versions)
      
       ch_SV_annotsv_vcf = SVPACK_ANNOTATE.out.tagged_vcf
   }

   //
   // ANNOT_SV for Somatic SV annotation
   //
   if ( !params.skip_severus && !params.skip_annotsv)
   {
        ch_SV_annotsv_cache = channel.of([[:],"${params.annotsv_cache}"])
        if ( !params.skip_annotsv_install )
        {
          ANNOTSV_INSTALLANNOTATIONS()
          ch_versions       = ch_versions.mix(ANNOTSV_INSTALLANNOTATIONS.out.versions)
          ch_SV_annotsv_cache = ANNOTSV_INSTALLANNOTATIONS.out.annotations
                                     .map { AnnotSV_annotations ->
                                      [[:], AnnotSV_annotations]
                                     }  
        }
        // ch_SV_annotsv_cache.view()      
     
        // TABIX_BGZIPTABIX(ch_SV_annotsv_vcf)
        TABIX_SV_VCF(ch_SV_annotsv_vcf)
        ch_versions       = ch_versions.mix(TABIX_SV_VCF.out.versions)
         
        ch_SV_annotsv_vcf_tbi =TABIX_SV_VCF.out.gz_tbi
                                .map { meta, vcf, tbi ->
                                  [meta, vcf, tbi, []] // [meta, sv_vcf, sv_vcf_tbi, candidate_small_variants]
                                }
        // ch_SV_annotsv_vcf_tbi.view()        
        
        ANNOTSV_ANNOTSV(ch_SV_annotsv_vcf_tbi, ch_SV_annotsv_cache, [[:],[]], [[:],[]], [[:],[]])
        ch_versions       = ch_versions.mix(ANNOTSV_ANNOTSV.out.versions)
   } 
     
   //
   //  CHORD using SNV and SV calling results
   //
   if( !params.skip_deepsomatic && !params.skip_severus && !params.skip_chord) {
        UNZIPFILES (ch_somatic_snv_vcf_gz)
        ch_chord_somatic_snv_vcf = UNZIPFILES.out.files
                             .map {pair_meta, snv_vcf ->
                             [pair_meta.id, pair_meta, snv_vcf]
                              }

        ch_chord_somatic_sv_vcf= ch_somatic_sv_vcf
                                 .map {pair_meta, sv_vcf ->
                                  [ pair_meta.id, pair_meta, sv_vcf ]
                                 } // [pair_id, pair_meta, somatic_vcf]        

        ch_chord_snv_sv_vcfs= ch_chord_somatic_snv_vcf.combine(ch_chord_somatic_sv_vcf, by :0)
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

        // ch_chord_snv_sv_vcfs.view()

        chord_genome_fasta  = ch_genome_fasta.map {meta, genome_fasta -> [genome_fasta] }
        chord_genome_fai    = ch_genome_fai.map { meta, genome_fai -> [genome_fai] }
        
        SAMTOOLS_DICT(ch_genome_fasta)
        ch_versions        = ch_versions.mix(SAMTOOLS_DICT.out.versions)
        
        chord_genome_dict =SAMTOOLS_DICT.out.dict
                           .map {meta, genome_dict -> [genome_dict] }

       // chord_genome_fasta.view()
       // chord_genome_fai.view()

        CHORD(ch_chord_snv_sv_vcfs, chord_genome_fasta, chord_genome_fai, chord_genome_dict)

        ch_versions        = ch_versions.mix(CHORD.out.versions)
   }
    
   //
   //  MUTATIONPATTERN for SNV mutatation signature 
   //
   if (!params.skip_mutationalpattern) {
       ch_mutationalpattern_vcf   = ch_somatic_snv_vcf_gz
       ch_mutationalpattern_genome= channel.of( [ [ id:'hg38' ], 'BSgenome.Hsapiens.UCSC.hg38' ] )   //  ch_genome_fasta
       MUTATIONALPATTERN(ch_mutationalpattern_vcf, ch_mutationalpattern_genome, params.mutationalpattern_max_delta)
       ch_versions        = ch_versions.mix(MUTATIONALPATTERN.out.versions)
   }
    
   //
   // CNVKIT for somatic CNV calling
   //
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
