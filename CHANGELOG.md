# nf-core/pacsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0 - 2026-01-24

Initial release of nf-core/pacsomatic, a comprehensive pipeline for somatic variant analysis from PacBio HiFi sequencing data.

### Features

- **Alignment**: Read alignment to reference genome using [pbmm2](https://github.com/PacificBiosciences/pbmm2)
- **Quality Control**: Alignment QC with SAMtools, mosdepth, and bamCoverage
- **Germline Variant Calling**: SNV/indel calling with [Clair3](https://github.com/HKU-BAL/Clair3)
- **Germline Phasing**: Haplotype phasing with [HiPhase](https://github.com/PacificBiosciences/HiPhase)
- **Somatic SNV/Indel Calling**: Deep learning-based calling with [DeepSomatic](https://github.com/google/deepsomatic)
- **Somatic Phasing**: Phase both germline and somatic variants in tumor samples with HiPhase
- **Variant Annotation**: Functional annotation with [Ensembl VEP](https://github.com/Ensembl/ensembl-vep)
- **Somatic SV Calling**: Structural variant detection with [Severus](https://github.com/KolmogorovLab/Severus)
- **SV Filtering/Annotation**: SV processing with [svpack](https://github.com/PacificBiosciences/svpack) and [AnnotSV](https://github.com/lgmgeo/AnnotSV)
- **CNV Calling**: Copy number variant detection with [CNVkit](https://github.com/etal/cnvkit)
- **Methylation Analysis**: CpG methylation calling with [pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools)
- **DMR Detection**: Differential methylation region detection with [DSS](https://bioconductor.org/packages/release/bioc/html/DSS.html)
- **DMR Annotation**: Annotation of DMRs with [annotatr](https://bioconductor.org/packages/release/bioc/html/annotatr.html)
- **Tumor Clonality**: Purity and ploidy estimation with [AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber), [COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt), and [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple)
- **Mutational Signatures**: Signature analysis with [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns)
- **HRD Estimation**: Homologous recombination deficiency prediction with [CHORD](https://github.com/UMCUGenetics/CHORD)
- **Reporting**: Aggregated QC reports with [MultiQC](http://multiqc.info/)

### Pipeline Features

- Flexible skip options for all major analysis steps
- Automatic BAM merging when multiple BAM files are provided per sample
- Tumor-normal paired analysis for all somatic calling
- Containerized execution with Docker, Singularity, and Conda support
- Full nf-core template compatibility

### Dependencies

- Nextflow >= 24.04.2
- nf-schema plugin 2.3.0
