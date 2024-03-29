# nf-core/crispriscreen: Usage

<!--
## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/crispriscreen/usage](https://nf-co.re/crispriscreen/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._
-->

## Introduction

The pipeline processes (Illumina) next generation sequencing data obtained from CRISPRi repression library screenings. It can also be used for other types of screenings as long as the input data is similar, that is, short PCR-amplified reads containing some sort of barcodes. The pipeline requires as little as three input types:

- Sequencing data in `*.fastq.gz` format
- The `samplesheet.csv` describing samples. See following section for details.
- A reference 'library' in `*.fasta` format

The results of the pipeline are described in detail in the [output](output.md) section.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A full samplesheet describing single-end sequencing data may look something like the one below. This is for 6 samples, 3 replicates each for two different time points of the same condition/treatment. In addition to a regular `samplesheet.csv`, the sheet for the `crispriscreen` pipeline takes the input columns `condition`, `replicate` , `time`, `group`, and `reference_group`. These columns are meta data describing the experimental design. This information is then used to calculate fitness scores which are a function of (generation) time. Each time point makes up one `group` and is compared to the `reference_group`, usually the zero time point of the same condition/treatment. If the additional columns are missing or incomplete, fitness score calculation is omitted but the pipeline will still be able to run. In this case the final output are read count tables from the `subread/featureCounts` module instead of fitness scores.

```console
sample,fastq_1,fastq_2,condition,replicate,time,group,reference_group
Example_S1_R1,./test/fastq/Example_S1_R1.fastq.gz,,example,1,0,1,1
Example_S1_R2,./test/fastq/Example_S1_R2.fastq.gz,,example,2,0,1,1
Example_S1_R3,./test/fastq/Example_S1_R3.fastq.gz,,example,3,0,1,1
Example_S2_R1,./test/fastq/Example_S2_R1.fastq.gz,,example,1,4,2,1
Example_S2_R2,./test/fastq/Example_S2_R2.fastq.gz,,example,2,4,2,1
Example_S2_R3,./test/fastq/Example_S2_R3.fastq.gz,,example,3,4,2,1
```

| Column            | Description                                                                                                                         |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `sample`          | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.                       |
| `fastq_1`         | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".          |
| `fastq_2`         | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".          |
| `condition`       | Name of the treatment or condition.                                                                                                 |
| `replicate`       | Index of the replicate.                                                                                                             |
| `date`            | Sampling date.                                                                                                                      |
| `time`            | Generation time or a similar metric after which library composition is sampled. Note: not meant to be a time stamp, but a duration! |
| `group`           | Number indicating all replicates of one treatment/time combination that belong together.                                            |
| `reference_group` | Number indicating the group that these samples will be compared to. For zero time points, can be the same as `group`.               |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```console
nextflow run <path/to/nf-core-crispriscreen> -profile <docker/singularity/podman/shifter/charliecloud/conda/institute> --input <sample_sheet> --fasta <fasta_file> --outdir <path/to/output>
```

To run the pipeline with the enclosed sample data in `assets/`, run:

```console
cd path/to/nf-core-crispriscreen
nextflow run ./ -profile singularity --input assets/samplesheet.csv" --fasta "assets/library.fasta" --outdir "results"
```

The pipeline was successfully tested with `docker` and `singularity` profiles. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Custom options to run the pipeline

The following command line options can be used to customize the behavior of the pipeline.

| Parameter                                     | Default value | Description                                                                                                                                                                                                                                                                                                                        |
| --------------------------------------------- | ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Subsampling** (`Seqtk/sample`)              |               |                                                                                                                                                                                                                                                                                                                                    |
| `subsampling`                                 | `false`       | Whether to subsample input reads.                                                                                                                                                                                                                                                                                                  |
| `sample_size`                                 | 10000         | Number of reads that are sampled per input file if subsampling is true.                                                                                                                                                                                                                                                            |
| `fixed_seed`                                  | `false`       | Whether the same reads (by position) should be sampled from input files. Can be useful for paired end reads.                                                                                                                                                                                                                       |
|                                               |               |                                                                                                                                                                                                                                                                                                                                    |
| **Generic quality trimming** (`Trim Galore!`) |               |                                                                                                                                                                                                                                                                                                                                    |
| `clip_r1`                                     | `null`        | Remove bp from the 5' end of read 1 (or single-end reads).                                                                                                                                                                                                                                                                         |
| `clip_r2`                                     | `null`        | Remove bp from the 5' end of read 2 (paired-end reads only).                                                                                                                                                                                                                                                                       |
| `three_prime_clip_r1`                         | `null`        | Remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed.                                                                                                                                                                                                                                             |
| `three_prime_clip_r2`                         | `null`        | Remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed.                                                                                                                                                                                                                                             |
| `trim_nextseq`                                | `null`        | Apply the --nextseq=X option, to trim based on quality after removing poly-G tails.                                                                                                                                                                                                                                                |
| `save_trimmed`                                | `false`       | Save the trimmed FastQ files in the results directory.                                                                                                                                                                                                                                                                             |
| `skip_trimming`                               | `false`       | Skip all adapter trimming (`Trim galore!` and `cutadapt`).                                                                                                                                                                                                                                                                         |
|                                               |               |                                                                                                                                                                                                                                                                                                                                    |
| **Specific read trimming** (`cutadapt`)       |               |                                                                                                                                                                                                                                                                                                                                    |
| `three_prime_adapter`                         | `null`        | The sequence of the 3' adapter / primer to be trimmed by cutadapt. This option has to be used for _linked_ adapters where both 3' and 5' adapters should be trimmed simultaneously (syntax: `--three_prime_adapter ADAPTER1...ADAPTER2`, see [cutadapt manual](https://cutadapt.readthedocs.io/en/stable/guide.html#basic-usage) ) |
| `five_prime_adapter`                          | `null`        | The sequence of the 5' adapter / primer to be trimmed by cutadapt.                                                                                                                                                                                                                                                                 |
| `error_rate`                                  | `false`       | Error rate for cutadapt. The default, 0.1, allows 1 mismatch per 10 nt.                                                                                                                                                                                                                                                            |
|                                               |               |                                                                                                                                                                                                                                                                                                                                    |
| **Alignment** (`Bowtie2`)                     |               |                                                                                                                                                                                                                                                                                                                                    |
| `save_unaligned`                              | `false`       | Whether to export unaligned reads from bowtie2 output or not.                                                                                                                                                                                                                                                                      |
| `sort_bam`                                    | `true`        | Whether to use samtools sort (true) or samtools view (false)                                                                                                                                                                                                                                                                       |
| `filter_mapq`                                 | `null`        | whether to filter reads by bowtie2's MAPQ value. The default is null (no filtering). (false)                                                                                                                                                                                                                                       |
|                                               |               |                                                                                                                                                                                                                                                                                                                                    |
| **Fitness calculation**                       |               |                                                                                                                                                                                                                                                                                                                                    |
| `normalization`                               | `false`       | If true, read counts for samples of identical time points are quantile normalized. Can be useful to counter act biases between different samples of the same time point, e.g. caused by under- or over-estimation of generation time.                                                                                              |
| `gene_fitness`                                | `true`        | Whether to calculate gene fitness from single sgRNA mutant or transposon fitness. Usually used for all types of barcode sequencing where several guide RNAs or transposons target the same gene.                                                                                                                                   |
| `gene_sep`                                    | `<pipe>`      | A string or symbol that separates target (gene) name from guide RNA index or transposon location in the fasta file. Usually used for all types of barcode sequencing where several guide RNAs or transposons target the same gene.                                                                                                 |
| `gene_controls`                               | `null`        | A string serving as regular expression to identify control barcodes/guide RNAs. It is matched against the target (gene) names in the fasta file. It is only used for p-value calculation using 2-sample Wilcoxon rank sum test. If omitted (the default), a 1-sample test is conducted                                             |
| `run_mageck`                                  | `true`        | Whether to run Mageck MLE for fitness calculation. For details, see [Mageck documentation](https://sourceforge.net/p/mageck/wiki/advanced_tutorial/#advanced-tutorials)                                                                                                                                                            |
| `design_matrix`                               | `null`        | Optional path to a design matrix file for Mageck fitness calculation. This parameter expects a tab-separated values file with mandatory columns `sample`, `baseline`, and one column each for each condition. For details, see [Mageck documentation](https://sourceforge.net/p/mageck/wiki/advanced_tutorial/#advanced-tutorials) |

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/crispriscreen
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/crispriscreen releases page](https://github.com/nf-core/crispriscreen/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/).

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.
- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/software/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
