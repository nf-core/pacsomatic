/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UTILITY FUNCTIONS FOR PACSOMATIC PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check required parameters are provided
//
def checkParameters() {
    def errors = []

    // Check required parameters
    if (!params.fasta) {
        errors << "Missing required parameter: --fasta (reference genome FASTA file)"
    }

    // Validate workflow type
    if (params.workflow && !(params.workflow in ['wgs', 'wes'])) {
        errors << "Invalid workflow type: '${params.workflow}'. Must be 'wgs' or 'wes'"
    }

    // VEP annotation requires all three parameters
    if (!params.skip_ensemblvep && !params.skip_deepsomatic) {
        if (!params.vep_assembly || !params.vep_species) {
            log.warn "VEP annotation requires --vep_assembly and --vep_species parameters. VEP will be skipped if not provided."
        }
    }

    // CNVkit WES mode requires target BED
    if (params.workflow == 'wes' && !params.skip_cnvkit && !params.cnv_target_bed) {
        log.warn "WES workflow with CNVkit requires --cnv_target_bed. CNVkit may fail without target regions."
    }

    // Tumor clonality analysis resource checks
    if (!params.skip_tumor_clonality) {
        def missingTumorClonalityResources = []
        if (!params.heterozygous_sites) missingTumorClonalityResources << 'heterozygous_sites'
        if (!params.gc_profile) missingTumorClonalityResources << 'gc_profile'
        if (!params.ensembl_data_dir) missingTumorClonalityResources << 'ensembl_data_dir'

        if (missingTumorClonalityResources) {
            log.warn "Tumor clonality analysis may require the following resources for optimal results: ${missingTumorClonalityResources.join(', ')}"
        }
    }

    // Methylation analysis requires phasing
    if (!params.skip_pbcpgtools && (params.skip_hiphase || params.skip_somatic_hiphase)) {
        log.warn "Methylation analysis requires both germline (--skip_hiphase=false) and somatic (--skip_somatic_hiphase=false) phasing to be enabled."
    }

    // Signature analysis requires both SNV and SV calling
    if ((!params.skip_chord || !params.skip_mutationalpattern) && (params.skip_deepsomatic || params.skip_sv)) {
        log.warn "Signature analysis (CHORD/MutationalPatterns) requires both DeepSomatic SNV calling and SV calling to be enabled."
    }

    // Exit with errors if any
    if (errors) {
        log.error "Parameter validation failed:"
        errors.each { error -> log.error "  - ${error}" }
        exit 1
    }
}

//
// Check path parameters exist
//
def checkPathParameters() {
    def pathParams = [
        [path: params.fasta, name: 'fasta'],
        [path: params.clair3_model_path, name: 'clair3_model_path'],
        [path: params.vep_cache, name: 'vep_cache'],
        [path: params.severus_trf_bed, name: 'severus_trf_bed'],
        [path: params.svpack_ctrl_vcf, name: 'svpack_ctrl_vcf'],
        [path: params.svpack_ref_gff, name: 'svpack_ref_gff'],
        [path: params.annotsv_cache, name: 'annotsv_cache'],
        [path: params.cnv_target_bed, name: 'cnv_target_bed'],
        [path: params.cnv_reference, name: 'cnv_reference'],
        [path: params.cnv_germline_vcf, name: 'cnv_germline_vcf'],
        [path: params.heterozygous_sites, name: 'heterozygous_sites'],
        [path: params.gc_profile, name: 'gc_profile'],
        [path: params.driver_gene_panel, name: 'driver_gene_panel'],
        [path: params.known_hotspots_somatic, name: 'known_hotspots_somatic'],
        [path: params.known_hotspots_germline, name: 'known_hotspots_germline'],
        [path: params.ensembl_data_dir, name: 'ensembl_data_dir'],
        [path: params.target_regions_bed, name: 'target_regions_bed'],
        [path: params.diploid_regions, name: 'diploid_regions'],
        [path: params.target_region_normalisation, name: 'target_region_normalisation']
    ]

    pathParams.each { param ->
        if (param.path) {
            try {
                file(param.path, checkIfExists: true)
            } catch (Exception e) {
                log.error "File not found for parameter --${param.name}: ${param.path}"
                exit 1
            }
        }
    }
}

//
// Validate input samplesheet structure
//
def validateSamplesheet(samplesheet) {
    def patients = [:]
    def errors = []

    samplesheet.each { row ->
        def patient = row.patient
        def sample = row.sample
        def status = row.status

        // Track samples per patient
        if (!patients.containsKey(patient)) {
            patients[patient] = [normal: [], tumor: []]
        }

        if (status == 0) {
            patients[patient].normal << sample
        } else if (status == 1) {
            patients[patient].tumor << sample
        }
    }

    // Validate each patient has at least one tumor and one normal
    patients.each { patient, samples ->
        if (samples.normal.size() == 0) {
            errors << "Patient '${patient}' has no normal sample (status=0)"
        }
        if (samples.tumor.size() == 0) {
            errors << "Patient '${patient}' has no tumor sample (status=1)"
        }
    }

    if (errors) {
        log.error "Samplesheet validation failed:"
        errors.each { error -> log.error "  - ${error}" }
        exit 1
    }

    return true
}
