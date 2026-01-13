<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-pacsomatic_logo_dark.png">
    <img alt="nf-core/pacsomatic" src="docs/images/nf-core-pacsomatic_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/pacsomatic/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/pacsomatic/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/pacsomatic/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/pacsomatic/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/pacsomatic/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/pacsomatic)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23pacsomatic-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/pacsomatic)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/pacsomatic** is a bioinformatics best-practice pipeline for somatic variant analysis using PacBio HiFi sequencing data from matched tumor-normal samples.

The pipeline performs comprehensive somatic analysis including:

- **Variant calling**: SNVs, indels, structural variants (SVs), and copy number variants (CNVs)
- **Variant annotation**: Functional annotation and mutation signatures
- **Methylation analysis**: CpG methylation calling and differential methylation region (DMR) detection
- **Tumor characterization**: Clonality, purity, ploidy, and homologous recombination deficiency (HRD) analysis

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow manager to run tasks across multiple compute infrastructures in a portable, reproducible manner. It is designed following the [nf-core](https://nf-co.re/) community's best practices and utilizes containerization with Docker, Singularity, or Conda for dependency management.

<p align="center">
    <img src="https://github.com/stjudecab/pacsomatic/blob/dev/docs/images/Pacsomatic_workflow_beta.png" alt="nf-core/pacsomatic workflow overview" width="80%"/>
</p>

## Pipeline Overview

The pipeline performs the following steps:

### 1. Quality Control and Alignment
- Read quality control ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](http://multiqc.info/))
- Read alignment to reference genome ([pbmm2](https://github.com/PacificBiosciences/pbmm2))
- Alignment sorting and indexing ([SAMtools](https://sourceforge.net/projects/samtools/files/samtools/))
- Alignment quality assessment ([bamCoverage](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html), [mosdepth](https://github.com/brentp/mosdepth))

### 2. Variant Calling
- **Germline SNVs**: [Clair3](https://github.com/HKU-BAL/Clair3)
- **Variant phasing**: [HiPhase](https://github.com/PacificBiosciences/HiPhase)
- **Somatic SNVs/indels**: [DeepSomatic](https://github.com/google/deepsomatic)
- **Somatic structural variants**: [Severus](https://github.com/KolmogorovLab/Severus)
- **Somatic copy number variants**: [CNVkit](https://github.com/etal/cnvkit)

### 3. Variant Annotation and Filtering
- SNV/indel functional annotation ([VEP](https://github.com/Ensembl/ensembl-vep))
- Mutation signature analysis ([MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns))
- SV filtering ([svpack](https://github.com/PacificBiosciences/svpack))
- SV annotation ([AnnotSV](https://github.com/lgmgeo/AnnotSV))

### 4. Methylation Analysis
- CpG methylation calling ([pb-CpG-tools](https://github.com/PacificBiosciences/pb-CpG-tools))
- Differential methylation region detection ([DSS](https://forge.irstea.fr/chloe.cerutti/bsseqmethdiffanalysis/-/blob/main/DSS/DMR.R))
- DMR annotation ([annotatr](https://bioconductor.org/packages/release/bioc/html/annotatr.html))

### 5. Tumor Characterization
- Homologous recombination deficiency estimation ([CHORD](https://github.com/UMCUGenetics/CHORD))
- Tumor purity and ploidy analysis ([AMBER](https://github.com/hartwigmedical/hmftools/tree/master/amber), [COBALT](https://github.com/hartwigmedical/hmftools/tree/master/cobalt), [PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

**Minimal samplesheet format** (`samplesheet.csv`):

```csv
patient,sample,status,bam
ID1,S1_tumor,1,/path/to/ID1_S1_tumor.bam
ID1,S1_normal,0,/path/to/ID1_S1_normal.bam
ID2,S2_tumor,1,/path/to/ID2_S2_tumor.bam
ID2,S2_normal,0,/path/to/ID2_S2_normal.bam
```

**Extended samplesheet with PBI index files**:

```csv
patient,sample,status,bam,pbi
ID1,S1_tumor,1,/path/to/ID1_S1_tumor.bam,/path/to/ID1_S1_tumor.bam.pbi
ID1,S1_normal,0,/path/to/ID1_S1_normal.bam,/path/to/ID1_S1_normal.bam.pbi
```

**Column descriptions**:
- `patient`: Unique patient identifier (samples with the same ID are treated as matched pairs)
- `sample`: Unique sample identifier
- `status`: Sample type (`1` = tumor, `0` = normal)
- `bam`: Full path to unaligned BAM file
- `pbi`: (Optional) Full path to PacBio index (.pbi) file

Now, you can run the pipeline using:

```bash
nextflow run nf-core/pacsomatic \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --genome GRCh38
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/pacsomatic/usage) and the [parameter documentation](https://nf-co.re/pacsomatic/parameters).

## Pipeline Output

Results are organized into functionally grouped subdirectories:

```
results/
├── alignment/                 # Aligned BAMs and QC metrics
│   ├── pbmm2/                # Aligned BAM files
│   └── qc/                   # Alignment quality control
├── germline_snv/             # Germline variants and phasing
│   ├── clair3/              # Germline SNV/indel calls
│   └── hiphase/             # Phased germline variants
├── somatic_snv/              # Somatic SNV/indel analysis
│   ├── deepsomatic/         # Somatic variant calls
│   ├── vep_annot/           # VEP annotations
│   └── hiphase_somatic/     # Phased somatic variants
├── somatic_sv/               # Structural variant analysis
│   ├── severus/             # SV calls
│   ├── svpack/              # Filtered SVs
│   └── annotsv_annot/       # SV annotations
├── somatic_cnv/              # Copy number variants
│   └── cnvkit/              # CNVkit results
├── methylation/              # Methylation analysis
│   ├── pb_cpg_tools/        # CpG methylation scores
│   ├── dss_dmr/             # Differential methylation regions
│   └── dmr_annot/           # DMR annotations
├── tumor_clonality/          # Tumor purity and ploidy
│   ├── amber/               # BAF analysis
│   ├── cobalt/              # Read depth ratios
│   └── purple/              # Purity/ploidy estimation
├── signature_analysis/       # Mutational signatures and HRD
│   ├── mutationalpattern/   # Mutation signatures
│   └── chord/               # HRD estimation
├── pipeline_info/            # Pipeline execution reports
└── multiqc/                  # Aggregated QC report
```

For detailed descriptions of output files, see the [output documentation](https://nf-co.re/pacsomatic/output).

To view example results from a full-size test dataset, visit the [results page](https://nf-co.re/pacsomatic/results) on the nf-core website.

## Credits

nf-core/pacsomatic was originally written by Wenchao Zhang and Haidong Yi.

We thank the following people for their extensive assistance in the development of this pipeline:

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#pacsomatic` channel](https://nfcore.slack.com/channels/pacsomatic) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/pacsomatic for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
