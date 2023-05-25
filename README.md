# ![nf-core/crispriscreen](docs/images/nf-core-crispriscreen_logo_light.png#gh-light-mode-only) ![nf-core/crispriscreen](docs/images/nf-core-crispriscreen_logo_dark.png#gh-dark-mode-only)

<!--
[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/crispriscreen/results)
[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
-->

[![GitHub Actions CI Status](https://github.com/MPUSP/nf-core-crispriscreen/workflows/nf-core%20CI/badge.svg)](https://github.com/MPUSP/nf-core-crispriscreen/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/MPUSP/nf-core-crispriscreen/workflows/nf-core%20linting/badge.svg?branch=dev)](https://github.com/MPUSP/nf-core-crispriscreen/actions?query=workflow%3A%22nf-core+linting%22)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/crispriscreen)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23crispriscreen-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/crispriscreen)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/crispriscreen** is a bioinformatics best-practice analysis pipeline to process next generation sequencing data obtained from CRISPRi repression library screenings.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable
On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/crispriscreen/results).
 -->

## Pipeline summary

1. Sub-sampling of reads ([`Seqtk/sample`](https://github.com/lh3/seqtk), optional)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Generic adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Specific primer sequence trimming ([`cutadapt`](https://cutadapt.readthedocs.io/en/stable/index.html))
5. Preparation of `*.fasta` library (custom [`R` script](https://cran.r-project.org/))
6. Alignment using ([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
   1. Build index from `*.fasta` library
   2. Align reads to library
   3. Optional filtering by mapping quality
7. Count reads per target and input file ([`subread/featurecounts`](https://nf-co.re/modules/subread_featurecounts))
8. Quantify gene fitness score from multiple targets per gene
   1. Option 1: Gene fitness is calculated using [`Mageck` MLE](https://sourceforge.net/p/mageck/wiki/Home/)
   2. Option 2: Gene fitness is calculated using [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
9. Generate HTML report with fitness results ([`R markdown`](https://nf-co.re/modules/rmarkdownnotebook))
10. Present QC for raw and mapped reads ([`MultiQC`](http://multiqc.info/))

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`).
   It is recommended to use [`Conda`](https://conda.io/miniconda.html) (or `mamba` / `micromamba`) to install all dependencies in a fresh environment.

   ```console
   conda create --name env_nf
   conda activate env_nf
   conda install -c conda-forge nextflow
   ```

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

   ```console
   conda install -c conda-forge singularity
   ```

3. Download the pipeline.

   ```console
   cd <your/target/dir>
   git clone https://github.com/MPUSP/nf-core-crispriscreen
   ```

4. Test it on the minimal dataset included with this repository. Since `nf-core-crispriscreen` is not a canonical `nf-core` pipeline (yet), it is necessary to indicate the path to the pipeline folder after the `run` statement.

   The generalized command to run the pipeline:

   ```console
   nextflow run <path/to/nf-core-crispriscreen> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input <sample_sheet> --fasta <fasta_file> --outdir <path/to/output>
   ```

   The command to run the pipeline on the enclosed test data using `Singularity` (recommended):

   ```console
   cd path/to/nf-core-crispriscreen
   nextflow run ./ -profile singularity --input "assets/samplesheet.csv" --fasta "assets/library.fasta" --outdir "results"
   ```

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

## Documentation

<!-- TODO: Add links to main nf-core website if published, e.g. [usage](https://nf-co.re/crispriscreen/usage) -->

The nf-core/crispriscreen pipeline comes with documentation about the pipeline [usage](https://MPUSP.github.io/nf-core-crispriscreen/usage) and [output](https://MPUSP.github.io/nf-core-crispriscreen/output).

## Credits

The following people contributed to the nf-core/crispriscreen pipeline:

- Dr. Michael Jahn, MPUSP Berlin ([ORCID](https://orcid.org/0000-0002-3913-153X)) - main author and maintainer
- Dr. Ute Hoffmann, Science for Life Lab, Stockholm - contributor

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#crispriscreen` channel](https://nfcore.slack.com/channels/crispriscreen) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

The following scientific articles are related to this pipeline.

The new version of CRISPRi libary for which this Nextflow pipeline was developed:

> Rui Miao, Michael Jahn, Kiyan Shabestary, Elton Paul Hudson.
> _CRISPR interference screens reveal tradeoffs between growth rate and robustness in Synechocystis sp. PCC 6803 across trophic conditions_
> bioRxiv 2023.02.13.528328. https://doi.org/10.1101/2023.02.13.528328

The original publication of the CRISPRi library for which a precursor of this pipeline was developed:

> Lun Yao, Kiyan Shabestary, Sarah Björk, Johannes Asplund-Samuelsson, Hakan Joensson, Michael Jahn & Elton Paul Hudson.
> _Pooled CRISPRi screening of the cyanobacterium Synechocystis sp PCC 6803 for enhanced industrial phenotypes_.
> Nature Communications, 11 (1666), 1–13. 2020. https://doi.org/10.1101/823534

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
