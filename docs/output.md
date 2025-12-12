# nf-core/pacsomatic: Output

## Introduction

This document describes the output produced by the Nextflow pipeline-nfcore/pacsomatic. Results will be stored in the directory specified during pipeline execution with the parameter: --outdir <OUTDIR>.

After the pipeline execution finishes, the output directory (<OUTDIR>) will contain several organized subdirectories, described below. All paths mentioned here are relative to the top-level output directory.
 
<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [PBMM2](#pbmm2) - Align samples to reference genome
- [Alignment QC](#alignment-qc) -  Calculate the QC metrics of each alignment  
   - [BAM_COVERAGE](#bamcoverage) - Generate the coverage tracks for each each alignment 
   - [MOSDEPTH](#mosdepth) - Calculate BAM depth for each alginemnt
   - [BAM_SORT_STATUS_SAMTOOLS](#bam_sort_status_samtools) -Use samtools sort/index/stats/flagstat for each alignment   
- [Methylation Detection and Annotation](#methylation-detection-and-annotation) Calculate CpG methylation for each alignment and detect and annotate the Differential Methylation Regions for tumor-normal pair
   - [CLAIR3](#clair3) - Germline SNV-INDEL  variant calling for all samples 
   - [HIPHASE](#hiphase) - Phase VCF and BAM files for normal samples
   - [HIPHASE_SOMATIC](#hiphase-somatic) - Phase VCF and BAM files for tumor samples
   - [PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES](#pbcpgtools-alignedbamtocpgscores) - Generate site methylation probabilities from mapped and phased BAM file
   - [DSS_DMR](#dss-dmr) - Use DSS(Dispersion Shrinkage for Sequencing) to detect DMR(Differential Methylation Region)    
   - [DMR_ANNOT](#dmr-annot) -Annotate the detected DMRs
- [Somatic CNV Calling](#somatic-cnv-calling) - Somatic CNV calling
   - [CNVKit](#cnvKit) - Use CNVKit packages to infer and visualize somatic copy number variants 
- [Somatic SNV INDEL Calling](#somatic-snv-indel-calling) - Somatic Variant call SNVs
   - [DEEPSOMATIC](#deepsomatic) -  Use deepsomatic to call somatic SNV and INDELs 
   - [VEP](#vep) - Use VEP for annotating somatic SNVs
   - [MUTATIONAL_PATTERN](#mutational_pattern) - Use Mutational_Patterns for mutation signature analysis
- [Somatic SV Calling](#somatic-sv-calling) - Somatic variant call SVs
   - [SEVERUS](#severus) - Use Severus to call somatic SVs
   - [SV_PACK](#sv-pack) - Use SV_Pack to filter called SVs
   - [ANNOT_SV](#annot-sv) - Use Annot-SV to annotate SVs    
- [Homologous Recombination Defficiency Estimation](#homologous-recombination-deficiency-estimation) -  Utilize the called SNVs and SVs to estimate the Homologous Recombination Defficiency (HRD).
   - [CHORD](#chord) - Use CHORD for HRD estimation 
- [Tumor purity and ploid estimation](#tumor-purity-and-ploid-estimation) - Estimate tumor purity and ploid
   - [AMBER](#amber) - Use hmftools -AMBER to analyze tumor-normal BAM pair to generate tumor BAF 
   - [COBALT](#cobalt) - Use htmlfools- COBALT to analyze tumor-normal BAM pair to determine the read depth ratio of the tumor against reference     
   - [PURPLE](#purple) - Use htmlfools- PURPLE to combine the BAF from AMBER and read depth ratio from COBALT to estimate tumor purity and ploid
- [FastQC](#fastqc) - Raw read QC- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution


## Output directory Structure
 ```
<OUTDIR>/
├── alignment_qc
├── genome
├── hrd_estimation
├── methylation_cpg
├── multiqc
├── pipeline_info
├── somatic_cnv
├── somatic_snv_indel
├── somatic_sv
└── tumor_purity_ploid
```

## PBMM2
<details markdown="1">
<summary>Output files</summary>

- `alignment_qc/pbmm2`
  - `<basename>.aligned.bam`: Aligned BAM

</details>
[PBMM2](https://github.com/PacificBiosciences/pbmm2) Aligned BAM files

## Alignment QC
### BAM_COVERAGE
<details markdown="1">
<summary>Output files</summary>

- `alignment_qc/bamcoverage`
  - `<basename>.bigWig`:  Covergare track

</details>
[bamcoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html) Generate coverage track

### MOSDEPTH
<details markdown="1">
<summary>Output files</summary>

- `alignment_qc/mosdepth`
  - `<basename>.mosdepth.global.dist.txt`: global cumulative distribution txt file indicating the proportion of total bases
  - `<basename>.mosdepth.summary.txt`:   mosdepth summary .txt file
  - `<basename>.per-base.bed.gz`:   mosdepth per-base .bed file
  - `<basename>.per-base.bed.gz.csi`:   mosdepth per-base .csi file

</details>
[mosdepth](https://github.com/brentp/mosdepth) Calculate depth 

### BAM_SORT_STATUS_SAMTOOLS
<details markdown="1">
<summary>Output files</summary>

- `alignment_qc/samtools`
  - `<basename>.bam`: samtool sort result file .bam
  - `<basename>.bam.bai`: samtool index result file .bam.bai
  - `<basename>.flagstat`: samtools flagstat result file .flagstat 
  - `<basename>.idxstats`: samtools idxstats result file .idxstats
  - `<basename>.stats`:   samtools stats result file .stats

</details>
[BAM_SORT_STATUS_SAMTOOLS](https://www.htslib.org/doc/samtools.html) use samtools to sort and index bams and generate some alignment QC metrics 

## Methylation Detection and Annotation
### clair3 
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/clair3`
  - `<basename>.vcf.gz`: germline variant file .vcf.gz 
  - `<basename>.vcf.gz.tbi`: germline variant tabix index file .vcf.gz.tbi

</details>
[clair3](https://github.com/HKU-BAL/Clair3) germline SNV calling for each alingment

### HIPHASE
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/hiphase`
  - `<basename>.phased.bam`: hiphase result file .phased.bam
  - `<basename>.phased.bam.bai`: hiphase result file .phased.bam.bai
  - `<basename>.phased.vcf`: hiphase result file .phased.vcf
  - `<basename>.stats.csv`: hiphase result file .stats.csv 

</details>
[HIPHASE](https://github.com/PacificBiosciences/HiPhase) Phasing the germline/normal alignment

### HIPHASE_SOMATIC
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/hiphase_somatic`
  - `<basename>.phased.bam`: hiphase_somatic result file .phased.bam
  - `<basename>.phased.bam.bai`: hiphase_somatic result file .phased.bam.bai          
  - `<basename>.germline_phased.vcf`: hiphas_somatice result file .germline_phased.vcf
  - `<basename>.somatic_phased.vcf`: hiphase_somatic result file .somatic_phased.vcf
  - `<basename>.stats.csv`: hiphase result file .stats.csv     

</details>
[HIPHASE_SOMATIC](https://github.com/PacificBiosciences/HiPhase) Phasing the tumor alignment

### PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES
This step apply to tumor and normal channel separatedly and deliver two corresponding sub-directories. 
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/pb_cpg_tools`
  - `<basename>.combined.bed.gz`: combined compressed bed file for complete read set
  - `<basename>.combined.bed.gz.tbi`: bed file index for complete read set  
  - `<basename>.combined.bw`: combined compressed .bw file for complete read set
  - `<basename>.hap1.bed.gz`: compressed bed file for haplotype 1
  - `<basename>.hap1.bed.gz.tbi`: bed file index for haplotype 1
  - `<basename>.hap1.bw`: compressed .bw file for haplotype 1
  - `<basename>.hap2.bed.gz`: compressed bed file for haplotype 2
  - `<basename>.hap2.bed.gz.tbi`: bed file index for haplotype 2
  - `<basename>.hap2.bw`: compressed .bw file for haplotype 2

</details>
[PBCPGTOOLS_ALIGNEDBAMTOCPGSCORES](https://github.com/PacificBiosciences/pb-CpG-tools) Pacbio CpG methylation tool to generate CpG score from aligned bam

### DSS_DMR
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/dss_dmr`
  - `<basename>.dmr.tsv`: differential methylation region .tsv file 

</details>
[DSS_DMR](https://rdrr.io/bioc/DSS/src/inst/doc/DSS.R) Use DSS (Dispersion shrinkage for sequencing) package to detect Differential Methylation Regions(DMR)  

### DMR_ANNOT
<details markdown="1">
<summary>Output files</summary>

- `methylation_cpg/dmr_annot`
  - `<prefix>_dmr_annotation_summary.tsv.gz`: DMR annotation-summary .tsv.gz file 
  - `<prefix>_hg38_genes_promoters_dmrs.tsv.gz`: DMR annotation- hg38 gene promoter .tsv.gz file 
  - `<prefix>_hg38_genes_1to5kb_dmrs.tsv.gz`: DMR annotation- hg38_genes_1to5kb .tsv.gz file
  - `<prefix>_hg38_genes_5UTRs_dmrs.tsv.gz`: DMR annotation- hg38_genes_5UTRs .tsv.gz file
  - `<prefix>_hg38_genes_exons_dmrs.tsv.gz`: DMR annotation-hg38_genes_exons 
  - `<prefix>_hg38_genes_introns_dmrs.tsv.gz`: DMR annotation- hg38_genes_introns
  - `<prefix>_hg38_genes_3UTRs_dmrs.tsv.gz`: DMR annotation- hg38_genes_3UTRs

</details>
[DMR_ANNOT](https://bioconductor.org/packages/release/bioc/html/annotatr.html) Annotate the Genomic Differential methylation Regions 

## Somatic CNV Calling
### CNVKit
<details markdown="1">
<summary>Output files</summary>

- `somatic_cnv/cnvkit`
  - `<prefix>_T.cnr`: .cnr file for tumor sample
  - `<prefix>_T.cns`: .cns file for tumor sample
  - `<prefix>_T-diagram.pdf`: .pdf file for tumor sample
  - `<prefix>_T-scatter.png`: .png file for tumor sample
  - `<prefix>_T_vs_<prefix>_N.cns`: Somatic CNV .cns file

</details>
[CNVKit](https://github.com/etal/cnvkit) CNVkit for CNV calling

## Somatic SNV INDEL Calling
### deepsomatic
<details markdown="1">
<summary>Output files</summary>

- `somatic_snv_indel/deepsomatic`
  - `<prefix>_T_vs_<prefix>_N.g.vcf.gz`: somatic SNV_INDEL variant .g.vcf.gz file
  - `<prefix>_T_vs_<prefix>_N.g.vcf.gz.tbi`: somatic SNV_INDEL variant .g.vcf.gz.tbi file
  - `<prefix>_T_vs_<prefix>_N.vcf.gz`: somatic SNV_INDEL variant .vcf.gz file
  - `<prefix>_T_vs_<prefix>_N.vcf.gz.tbi`: somatic SNV_INDEL variant .vcf.gz.tbi file

</details>
[DEEPSOMATIC](https://github.com/google/deepsomatic) Deepsomatic for somatic SNV_INDEL calling

### VEP
<details markdown="1">
<summary>Output files</summary>

- `somatic_snv_indel/vep_annot`
  - `<prefix>_T_vs_<prefix>_N.vep.anno.vcf.gz`: VEP annotated somatic SNV_INDEL variant .vep.anno.vcf.gz file
  - `<prefix>_T_vs_<prefix>_N.vep.anno.vcf.gz.tbi`: VEP annotated somatic SNV_INDEL variant .vep.anno.vcf.gz.tbi file
  - `<prefix>_T_vs_<prefix>_N.vep.anno.vcf.gz_summary.html`: VEP annotation summary .vep.anno.vcf.gz_summary.html 

</details>
[VEP](https://github.com/Ensembl/ensembl-vep) VEP annotation for somatic SNV_INDEL variants

### mutational_pattern
<details markdown="1">
<summary>Output files</summary>

- `somatic_snv_indel/mutationalpattern`
  - `<prefix>_T_vs_<prefix>_N.mutation_profile.pdf`: mutationalpattern .pdf file
  - `<prefix>_T_vs_<prefix>_N.mut_sigs_bootstrapped.tsv`: mutational signature boosttrapped .tsv file
  - `<prefix>_T_vs_<prefix>_N.mut_sigs.tsv`: mutational signature .tsv file
  - `<prefix>_T_vs_<prefix>_N.reconstructed_sigs.tsv`: mutational reconstructed signature .tsv file
  - `<prefix>_T_vs_<prefix>_N.type_occurences.tsv`: mutational signature type occurences .tsv file

</details>
[MUTATIONAL_PATTERN](https://github.com/UMCUGenetics/MutationalPatterns) Mutational Pattern annotation for mutational signature analysis

## Somatic SV Calling
### Severus
<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/severus/<prefix>`
  - `severus.log`: Severus run log file
  - `read_qual.txt`: Read quality .txt file
  - `breakpoints_double.csv`: breakpoint double .csv file
  - `all_SVs/severus_all.vcf`: All SVs called by Severus
  - `somatic_SVs/severus_somatic.vcf`: Somatic SVs called by Severus

</details>
[SEVERUS](https://github.com/KolmogorovLab/Severus) Somatic SV calling

### SV_Pack
<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/SVPack/<prefix>`
  - `SVPACK_FILTER.out.vcf`: SVPack filter result .out.vcf file
  - `SVPACK_MATCH.out.vcf`:  SVPack match result .out.vcf file 
  - `SVPACK_CONSEQUENCE.out.vcf`: SVPack consequence result .out.vcf file 
  - `SVPACK_TAGZYGOSITY.out.vcf`: SVPack tagzygosity result .out.vcf file

</details>
[SV_PACK](https://github.com/PacificBiosciences/svpack) filtering, comparing, and annotating SV

### ANNOT_SV
<details markdown="1">
<summary>Output files</summary>

- `somatic_sv/annotsv_annot`
  - `<basename>.tsv`: annot_sv annotation result .tsv file

</details>
[SV_PACK](https://github.com/lgmgeo/AnnotSV) annotation of SV

## Homologous Recombination Deficiency Estimation
### CHORD
<details markdown="1">
<summary>Output files</summary>

- `hrd_estimation/chord`
  - `<basename>.chord.mutation_contexts.tsv`: chord mutation context .tsv file
  - `<basename>.chord.prediction.tsv`: chord prediction result .tsv file

</details>
[CHORD](https://github.com/UMCUGenetics/CHORD) use chord for HRD Estimation

## Tumor purity and ploid estimation
### AMBER
<details markdown="1">
<summary>Output files</summary>

- `tumor_purity_ploid/amber`
  - `amber.version`: amber version
  - `<basename>.amber.homozygousregion.tsv`:  amber result .tsv file
  - `<basename>.amber.baf.pcf`: amber result .baf.pcf file
  - `<basename>.amber.baf.tsv.gz`: amber result .baf.tsv.gz file
  - `<basename>.amber.contamination.tsv`: amber result .contamination.tsv file
  - `<basename>.amber.contamination.vcf.gz`: amber result .contamination.vcf.gz file
  - `<basename>.amber.contamination.vcf.gz.tbi`: amber result .contamination.vcf.gz.tbi file 
  - `<basename>.amber.qc`: amber QC result .qc file 

</details>
[AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber) use hmftools-amber to generate tumor BAF file

### COBALT
<details markdown="1">
<summary>Output files</summary>

- `tumor_purity_ploid/cobalt`
  - `cobalt.version`: cobalt version file
  - `<basename>.cobalt.gc.median.tsv`: cobalt result .gc.median.tsv file
  - `<basename>.cobalt.ratio.median.tsv`:  cobalt result .ratio.median .tsv file
  - `<basename>.cobalt.ratio.pcf`: cobalt result .ratio.pcf file
  - `<basename>.cobalt.gc.median.tsv`: cobalt result .gc.median.tsv file

</details>
[COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt) use hmftools-cobalt to determine read depth ratios of tumor and reference/normal genomes

### PURPLE
<details markdown="1">
<summary>Output files</summary>

- `tumor_purity_ploid/purple`
  - `purple.version`: purple version file
  - `<prefix>_T.purple.cnv.gene.tsv`: purple result .cnv.gene.tsv file
  - `<prefix>_T.purple.cnv.somatic.tsv`: purple result .cnv.somatic.tsv file
  - `<PREFIX>_T.purple.driver.catalog.germline.tsv`: purple result .driver.catalog.germline.tsv file
  - `<PREFIX>_T.purple.driver.catalog.somatic.tsv`:  purple result .driver.catalog.somatic.tsv file
  - `<PREFIX>_T.purple.germline.deletion.tsv`: purple result .germline.deletion.tsv file
  - `<PREFIX>_T.purple.purity.range.tsv`: purple result .purity.range.tsv file
  - `<PREFIX>_T.purple.purity.tsv`: purple result tumor purity.tsv file
  - `<PREFIX>_T.purple.qc`: purple result QC file
  - `<PREFIX>_T.purple.segment.tsv`: purple result .segment.tsv file
  - `<PREFIX>_T.purple.somatic.clonality.tsv`: purple result .somatic.clonality.tsv file
  - `plot/<PREFIX>_T.purity.range.png `: purple result .purity.range.png file
  - `plot/<PREFIX>_T.segment.png `: purple result .segment.png file
  - `plot/<PREFIX>_T.somatic_data.tsv`: purple result .somatic_data.tsv file

</details>
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) use hmftools-purple to estimate purity ploidy

### FastQC
### 

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
