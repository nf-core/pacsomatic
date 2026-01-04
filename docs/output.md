# nf-core/pacsomatic: Output

## Introduction

This document describes the output produced by the Nextflow pipeline-nfcore/pacsomatic. Results will be stored in the directory specified during pipeline execution with the parameter: --outdir <OUTDIR>.

After the pipeline execution finishes, the output directory (<OUTDIR>) will contain several organized subdirectories, described below. All paths mentioned here are relative to the top-level output directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [PBMM2](#pbmm2) - Align samples to reference genome
- [Alignment QC](#alignment-qc) - Calculate the QC metrics of each alignment
   - [BAM_COVERAGE](#bam_coverage) - Generate the coverage tracks for each alignment
   - [MOSDEPTH](#mosdepth) - Calculate BAM depth for each alignment
   - [BAM_SORT_STATUS_SAMTOOLS](#bam_sort_status_samtools) - Use samtools sort/index/stats/flagstat for each alignment
- [Methylation Detection and Annotation](#methylation-detection-and-annotation) - Calculate CpG methylation for each alignment and detect and annotate the Differential Methylation Regions for tumor-normal pair
   - [Clair3](#clair3) - Germline SNV-INDEL variant calling for all samples
   - [HiPhase](#hiphase) - Phase VCF and BAM files for normal samples
   - [HiPhase Somatic](#hiphase_somatic) - Phase VCF and BAM files for tumor samples
   - [PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES](#pbcpgtools_alignedbamtocpgscores) - Generate site methylation probabilities from mapped and phased BAM file
   - [DSS_DMR](#dss_dmr) - Use DSS (Dispersion Shrinkage for Sequencing) to detect DMR (Differential Methylation Region)
   - [DMR_ANNOT](#dmr_annot) - Annotate the detected DMRs
- [Somatic CNV Calling](#somatic-cnv-calling) - Somatic CNV calling
   - [CNVkit](#cnvkit) - Use CNVkit packages to infer and visualize somatic copy number variants
- [Somatic SNV INDEL Calling](#somatic-snv-indel-calling) - Somatic Variant call SNVs
   - [DeepSomatic](#deepsomatic) - Use deepsomatic to call somatic SNV and INDELs
   - [VEP](#vep) - Use VEP for annotating somatic SNVs
   - [Mutational Pattern](#mutational_pattern) - Use Mutational_Patterns for mutation signature analysis
- [Somatic SV Calling](#somatic-sv-calling) - Somatic variant call SVs
   - [Severus](#severus) - Use Severus to call somatic SVs
   - [SV_Pack](#sv_pack) - Use SV_Pack to filter called SVs
   - [AnnotSV](#annot_sv) - Use AnnotSV to annotate SVs
- [Homologous Recombination Deficiency Estimation](#homologous-recombination-deficiency-estimation) - Utilize the called SNVs and SVs to estimate the Homologous Recombination Deficiency (HRD).
   - [CHORD](#chord) - Use CHORD for HRD estimation
- [Tumor purity and ploid estimation](#tumor-purity-and-ploid-estimation) - Estimate tumor purity and ploid
   - [AMBER](#amber) - Use hmftools-AMBER to analyze tumor-normal BAM pair to generate tumor BAF
   - [COBALT](#cobalt) - Use hmftools-COBALT to analyze tumor-normal BAM pair to determine the read depth ratio of the tumor against reference
   - [PURPLE](#purple) - Use hmftools-PURPLE to combine the BAF from AMBER and read depth ratio from COBALT to estimate tumor purity and ploid
- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution


## Output directory Structure

```
<OUTDIR>/
├── alignment/              # Aligned BAMs and QC metrics
│   ├── pbmm2/             # Aligned and sorted BAM files
│   └── qc/                # Alignment quality control
│       ├── samtools/      # SAMtools statistics
│       ├── bamcoverage/   # Coverage tracks (bigWig files)
│       └── mosdepth/      # Depth coverage analysis
├── genome/                 # Reference genome files
├── germline_snv/          # Germline variant calling and phasing
│   ├── clair3/           # Germline SNV/indel calls
│   └── hiphase/          # Phased germline variants and BAMs
├── somatic_snv/           # Somatic SNV/indel analysis
│   ├── deepsomatic/      # Somatic variant calls
│   ├── vep_annot/        # VEP annotations
│   └── hiphase_somatic/  # Phased somatic variants
├── somatic_sv/            # Somatic structural variant analysis
│   ├── severus/          # SV calls from Severus
│   ├── svpack/           # Filtered and annotated SVs
│   └── annotsv_annot/    # AnnotSV annotations
├── somatic_cnv/           # Somatic copy number variants
│   └── cnvkit/           # CNVkit results
├── methylation/           # Methylation analysis
│   ├── pb_cpg_tools/     # CpG methylation scores
│   ├── dss_dmr/          # Differential methylation regions
│   └── dmr_annot/        # DMR annotations
├── tumor_clonality/       # Tumor purity and ploidy
│   ├── amber/            # BAF analysis
│   ├── cobalt/           # Read depth ratios
│   └── purple/           # Purity and ploidy estimation
├── signature_analysis/    # Mutational signatures and HRD
│   ├── mutationalpattern/# Mutation signatures
│   └── chord/            # HRD estimation
├── multiqc/               # Aggregated QC report
└── pipeline_info/         # Pipeline execution reports
```

## PBMM2

<details markdown="1">
<summary>Output files</summary>

- `alignment/pbmm2/<sample_id>/`
  - `<sample_id>.aligned.bam`: Aligned and sorted BAM file
  - `<sample_id>.aligned.bam.bai`: BAM index file

</details>

[PBMM2](https://github.com/PacificBiosciences/pbmm2) aligns PacBio HiFi reads to the reference genome using minimap2.

## Alignment QC

### BAM_COVERAGE
<a id="bam_coverage"></a>

<details markdown="1">
<summary>Output files</summary>

- `alignment/qc/bamcoverage/<sample_id>/`
  - `<sample_id>.bigWig`: Coverage track in bigWig format

</details>

[bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) generates coverage tracks for visualization in genome browsers.

### MOSDEPTH

<details markdown="1">
<summary>Output files</summary>

- `alignment/qc/mosdepth/<sample_id>/`
  - `<sample_id>.mosdepth.global.dist.txt`: Global cumulative distribution
  - `<sample_id>.mosdepth.summary.txt`: Depth summary statistics
  - `<sample_id>.per-base.bed.gz`: Per-base coverage
  - `<sample_id>.per-base.bed.gz.csi`: Index file

</details>

[mosdepth](https://github.com/brentp/mosdepth) calculates depth coverage statistics for each alignment.

### BAM_SORT_STATUS_SAMTOOLS

<details markdown="1">
<summary>Output files</summary>

- `alignment/qc/samtools/<sample_id>/`
  - `<sample_id>.flagstat`: Alignment flag statistics
  - `<sample_id>.idxstats`: Index statistics
  - `<sample_id>.stats`: Detailed alignment statistics

</details>

[SAMtools](https://www.htslib.org/doc/samtools.html) generates comprehensive alignment quality metrics.

## Methylation Detection and Annotation

### Clair3

<details markdown="1">
<summary>Output files</summary>

- `germline_snv/clair3/<sample_id>/`
  - `<sample_id>_pileup.vcf.gz`: Germline variant calls
  - `<sample_id>_pileup.vcf.gz.tbi`: Tabix index

</details>

[Clair3](https://github.com/HKU-BAL/Clair3) performs germline SNV/indel calling for each sample.

### HiPhase

<details markdown="1">
<summary>Output files</summary>

- `germline_snv/hiphase/<sample_id>/`
  - `<sample_id>.phased.bam`: Phased BAM file
  - `<sample_id>.phased.bam.bai`: BAM index
  - `<sample_id>.phased.vcf.gz`: Phased variants
  - `<sample_id>.stats.csv`: Phasing statistics

</details>

[HiPhase](https://github.com/PacificBiosciences/HiPhase) phases germline variants and reads for normal samples.

### HiPhase Somatic
<a id="hiphase_somatic"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_snv/hiphase_somatic/<sample_id>/`
  - `<sample_id>.phased.bam`: Phased tumor BAM
  - `<sample_id>.phased.bam.bai`: BAM index
  - `<sample_id>.germline_phased.vcf.gz`: Phased germline variants
  - `<sample_id>.somatic_phased.vcf.gz`: Phased somatic variants
  - `<sample_id>.stats.csv`: Phasing statistics

</details>

[HiPhase](https://github.com/PacificBiosciences/HiPhase) phases both germline and somatic variants in tumor samples.

### PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES
<a id="pbcpgtools_alignedbamtocpgscores"></a>

This analysis is performed separately for tumor and normal samples.

<details markdown="1">
<summary>Output files</summary>

- `methylation/pb_cpg_tools/normal/<sample_id>/`
- `methylation/pb_cpg_tools/tumor/<sample_id>/`
  - `<sample_id>.combined.bed.gz`: Combined CpG scores (all reads)
  - `<sample_id>.combined.bed.gz.tbi`: Index file
  - `<sample_id>.combined.bw`: Coverage track (all reads)
  - `<sample_id>.hap1.bed.gz`: CpG scores for haplotype 1
  - `<sample_id>.hap1.bed.gz.tbi`: Index file
  - `<sample_id>.hap1.bw`: Coverage track for haplotype 1
  - `<sample_id>.hap2.bed.gz`: CpG scores for haplotype 2
  - `<sample_id>.hap2.bed.gz.tbi`: Index file
  - `<sample_id>.hap2.bw`: Coverage track for haplotype 2

</details>

[pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools) generates CpG methylation scores from phased BAM files.

### DSS_DMR
<a id="dss_dmr"></a>

<details markdown="1">
<summary>Output files</summary>

- `methylation/dss_dmr/<patient_id>/`
  - `<tumor_vs_normal>.dmr.tsv`: Differential methylation regions

</details>

[DSS](https://rdrr.io/bioc/DSS/src/inst/doc/DSS.R) detects differential methylation regions (DMRs) between tumor and normal samples.

### DMR_ANNOT
<a id="dmr_annot"></a>

<details markdown="1">
<summary>Output files</summary>

- `methylation/dmr_annot/<patient_id>/`
  - `<patient_id>_dmr_annotation_summary.tsv.gz`: DMR annotation summary
  - `<patient_id>_hg38_genes_promoters_dmrs.tsv.gz`: DMRs in gene promoters
  - `<patient_id>_hg38_genes_1to5kb_dmrs.tsv.gz`: DMRs 1-5kb from genes
  - `<patient_id>_hg38_genes_5UTRs_dmrs.tsv.gz`: DMRs in 5'UTRs
  - `<patient_id>_hg38_genes_exons_dmrs.tsv.gz`: DMRs in exons
  - `<patient_id>_hg38_genes_introns_dmrs.tsv.gz`: DMRs in introns
  - `<patient_id>_hg38_genes_3UTRs_dmrs.tsv.gz`: DMRs in 3'UTRs

</details>

[annotatr](https://bioconductor.org/packages/release/bioc/html/annotatr.html) annotates DMRs with genomic features.

## Somatic CNV Calling

### CNVkit
<a id="cnvkit"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_cnv/cnvkit/batch/<patient_id>/`
  - `<tumor_id>.cnr`: Copy number ratios
  - `<tumor_id>.cns`: Copy number segments
  - `<tumor_id>-diagram.pdf`: Chromosome diagram
  - `<tumor_id>-scatter.png`: Scatter plot
- `somatic_cnv/cnvkit/call/<patient_id>/`
  - `<tumor_id>.call.cns`: Called copy number segments

</details>

[CNVkit](https://github.com/etal/cnvkit) infers and visualizes somatic copy number variants from tumor-normal pairs.

## Somatic SNV INDEL Calling

### DeepSomatic
<a id="deepsomatic"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_snv/deepsomatic/<patient_id>/`
  - `<tumor_vs_normal>.g.vcf.gz`: GVCF with all sites
  - `<tumor_vs_normal>.g.vcf.gz.tbi`: Tabix index
  - `<tumor_vs_normal>.vcf.gz`: Somatic variants
  - `<tumor_vs_normal>.vcf.gz.tbi`: Tabix index

</details>

[DeepSomatic](https://github.com/google/deepsomatic) calls somatic SNVs and indels using deep learning.

### VEP

<details markdown="1">
<summary>Output files</summary>

- `somatic_snv/vep_annot/`
  - `<tumor_vs_normal>.vep.anno.vcf.gz`: VEP annotated variants
  - `<tumor_vs_normal>.vep.anno.vcf.gz.tbi`: Tabix index
  - `<tumor_vs_normal>.vep.anno.vcf.gz_summary.html`: Annotation summary

</details>

[VEP](https://github.com/Ensembl/ensembl-vep) provides functional annotation for somatic variants.

### Mutational Pattern
<a id="mutational_pattern"></a>

<details markdown="1">
<summary>Output files</summary>

- `signature_analysis/mutationalpattern/<patient_id>/`
  - `<tumor_vs_normal>.mutation_profile.pdf`: Mutation profile plots
  - `<tumor_vs_normal>.mut_sigs_bootstrapped.tsv`: Bootstrapped signatures
  - `<tumor_vs_normal>.mut_sigs.tsv`: Mutational signatures
  - `<tumor_vs_normal>.reconstructed_sigs.tsv`: Reconstructed signatures
  - `<tumor_vs_normal>.type_occurences.tsv`: Mutation type occurrences

</details>

[MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) identifies mutational signatures from somatic variants.

## Somatic SV Calling

### Severus
<a id="severus"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/severus/<patient_id>/orig/`
  - `<patient_id>_severus_somatic_SVs/severus_somatic.vcf`: Somatic SVs
  - `<patient_id>_severus_all_SVs/severus_all.vcf`: All SVs (somatic + germline)
  - `severus.log`: Run log
  - `read_qual.txt`: Read quality metrics
  - `breakpoints_double.csv`: Breakpoint details
- `somatic_sv/severus/<patient_id>/filtered/`
  - `<patient_id>.severus_somatic.vcf.gz`: Filtered somatic SVs (if svpack enabled)
  - `<patient_id>.severus_somatic.vcf.gz.tbi`: Tabix index

</details>

[Severus](https://github.com/KolmogorovLab/Severus) identifies somatic structural variants from tumor-normal pairs.

### SV_Pack
<a id="sv_pack"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/svpack/<patient_id>/`
  - `SVPACK_FILTER.out.vcf`: Filtered SVs
  - `SVPACK_MATCH.out.vcf`: Matched against control panel
  - `SVPACK_CONSEQUENCE.out.vcf`: Functional consequences
  - `SVPACK_TAGZYGOSITY.out.vcf`: Zygosity annotations

</details>

[svpack](https://github.com/PacificBiosciences/svpack) filters and annotates structural variants.

### AnnotSV
<a id="annot_sv"></a>

<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/annotsv_annot/<patient_id>/`
  - `<patient_id>.tsv`: Comprehensive SV annotations

</details>

[AnnotSV](https://github.com/lgmgeo/AnnotSV) provides detailed annotations for structural variants.

## Homologous Recombination Deficiency Estimation

### CHORD
<a id="chord"></a>

<details markdown="1">
<summary>Output files</summary>

- `signature_analysis/chord/<patient_id>/`
  - `<tumor_vs_normal>.chord.mutation_contexts.tsv`: Mutation contexts
  - `<tumor_vs_normal>.chord.prediction.tsv`: HRD predictions

</details>

[CHORD](https://github.com/UMCUGenetics/CHORD) estimates homologous recombination deficiency from mutational signatures.

## Tumor purity and ploid estimation

### AMBER
<a id="amber"></a>

<details markdown="1">
<summary>Output files</summary>

- `tumor_clonality/<patient_id>/amber/`
  - `amber.version`: Tool version
  - `<tumor_id>.amber.homozygousregion.tsv`: Homozygous regions
  - `<tumor_id>.amber.baf.pcf`: BAF piecewise constant fit
  - `<tumor_id>.amber.baf.tsv.gz`: B-allele frequencies
  - `<tumor_id>.amber.contamination.tsv`: Contamination estimates
  - `<tumor_id>.amber.contamination.vcf.gz`: Contamination VCF
  - `<tumor_id>.amber.contamination.vcf.gz.tbi`: Index file
  - `<tumor_id>.amber.qc`: Quality control metrics

</details>

[AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) calculates B-allele frequencies for tumor purity estimation.

### COBALT
<a id="cobalt"></a>

<details markdown="1">
<summary>Output files</summary>

- `tumor_clonality/<patient_id>/cobalt/`
  - `cobalt.version`: Tool version
  - `<tumor_id>.cobalt.gc.median.tsv`: GC-corrected median ratios
  - `<tumor_id>.cobalt.ratio.median.tsv`: Read depth ratios
  - `<tumor_id>.cobalt.ratio.pcf`: Piecewise constant fit
  - `<tumor_id>.cobalt.chr.len.tsv`: Chromosome lengths

</details>

[COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) calculates read depth ratios between tumor and normal.

### PURPLE
<a id="purple"></a>

<details markdown="1">
<summary>Output files</summary>

- `tumor_clonality/<patient_id>/purple/`
  - `purple.version`: Tool version
  - `<tumor_id>.purple.cnv.gene.tsv`: Gene-level CNV calls
  - `<tumor_id>.purple.cnv.somatic.tsv`: Somatic CNV segments
  - `<tumor_id>.purple.driver.catalog.germline.tsv`: Germline drivers
  - `<tumor_id>.purple.driver.catalog.somatic.tsv`: Somatic drivers
  - `<tumor_id>.purple.germline.deletion.tsv`: Germline deletions
  - `<tumor_id>.purple.purity.range.tsv`: Purity solution range
  - `<tumor_id>.purple.purity.tsv`: Final purity and ploidy estimates
  - `<tumor_id>.purple.qc`: Quality control metrics
  - `<tumor_id>.purple.segment.tsv`: Copy number segments
  - `<tumor_id>.purple.somatic.clonality.tsv`: Variant clonality
  - `plot/<tumor_id>.purity.range.png`: Purity range plot
  - `plot/<tumor_id>.segment.png`: Copy number plot
  - `plot/<tumor_id>.somatic_data.tsv`: Plot data

</details>

[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) estimates tumor purity and ploidy by combining AMBER and COBALT results.

### FastQC
<a id="fastqc"></a>



<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MultiQC
<a id="multiqc"></a>

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information
<a id="pipeline-information"></a>

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
