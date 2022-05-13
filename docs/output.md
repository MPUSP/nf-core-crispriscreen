# nf-core/crispriscreen: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

1. Sub-sampling of reads ([`Seqtk/sample`](https://github.com/lh3/seqtk), optional)
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
4. Preparation of `*.fasta` library (custom [`R` script](https://cran.r-project.org/))
5. Alignment using ([`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
   1. Build index from `*.fasta` library
   2. Align reads to library
6. Count reads per target and input file ([`subread/featurecounts`](https://nf-co.re/modules/subread_featurecounts))
7. Quantify gene fitness score from multiple targets per gene, report statistics ([`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
8. Generate HTML report with fitness results ([`R markdown`](https://nf-co.re/modules/rmarkdownnotebook))
9. Present QC for raw and mapped reads ([`MultiQC`](http://multiqc.info/))

### Seqtk/Sample

<details markdown="1">
<summary>Output files</summary>

- `seqtk/`
  - `*.subsampled.fastq.gz`: Subsampled compressed `fastq.gz` sequencing reads.

</details>

Sub-sampling of reads ([Seqtk/sample](https://github.com/lh3/seqtk), optional).

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

![MultiQC - FastQC sequence counts plot](images/mqc_fastqc_counts.png)

![MultiQC - FastQC mean quality scores plot](images/mqc_fastqc_quality.png)

![MultiQC - FastQC adapter content plot](images/mqc_fastqc_adapter.png)

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality.

### Trimgalore

<details markdown="1">
<summary>Output files</summary>

- `trimgalore/`
  - `*_trimming_report.txt`: Report of the trimming results for each `*.fastq.gz` input file.

</details>

Adapter and quality trimming of reads ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)).

### Preparation of `*.fasta` library

<details markdown="1">
<summary>Output files</summary>

- `prepare/`
  - `*.saf`: The input library in `saf` format converted from the provided `.fasta` file.
  - `*_controls.tsv`: _Optional_ table in tab-separated format with overview of control barcodes.

</details>

This module generates a `.tsv` input table from the provided `.fasta` file.
It also checks if a pattern for control barcodes has been passed with the `--gene_controls` parameter.
In case a pattern has been supplied but no matching barcodes were found, it stops with an error.

### Bowtie2

<details markdown="1">
<summary>Output files</summary>

- `bowtie2/bowtie2/`
  - `*.bt2`: Bowtie2 index created from the libary `*.fasta` file.
- `bowtie2/`
  - `*.bam`: Compressed sequence alignment files, one per input `.fastq.gz`.
  - `*.bowtie2.log`: Bowtie2 log file, one per input `.fastq.gz`.
  - `*.unmapped.fastq.gz`: Optionally exported unmapped reads, one per input `.fastq.gz`.

</details>

[`Bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) is used for mapping reads to the 'genome', here the library of guide RNAs/barcodes.

### FeatureCounts

<details markdown="1">
<summary>Output files</summary>

- `subread/`
  - `*.txt`: Main result of this module, a table with detailed read counts per target (guide RNA/barcode)
  - `*.txt.summary`: Summary of mapped and unmapped reads.

</details>

Summarizes read counts per target and input file, see [`subread/featurecounts`](https://nf-co.re/modules/subread_featurecounts).
This is the input for fitness score calculation with DESeq2.

### DESeq2

<details markdown="1">
<summary>Output files</summary>

- `fitness/`
  - `all_counts.tsv`: A summary table with all read counts per target (gene, barcode, sgRNA, ...), concatenated from the individual [#featurecounts] output files.
  - `result.tsv`: Table with fitness scores and other statistics for all conditions in `*.tsv` format.
  - `result.Rdata`: Table with fitness scores and other statistics for all conditions in memory-efficient `*.Rdata` format. Can be read into `R` using `load("result.Rdata`).

</details>

A custom R script employing [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) to quantify gene fitness score from multiple targets per gene, reporting different summary statistics. The final output of this module is a table in `*.txt` and `*.Rdata` format with the following columns:

| Column               | Type      | Example      | Comment                                                 |
| -------------------- | --------- | ------------ | ------------------------------------------------------- |
| sgRNA                | `chr`     | `aat_111`    | name of sgRNA as in `.fasta` reference                  |
| sgRNA_target         | `chr`     | `aat`        | name of sgRNA target gene                               |
| sgRNA_position       | `numeric` | 111          | position of sgRNA relative to target start              |
| condition            | `chr`     | `example`    | experimental condition                                  |
| date                 | `chr`     | `2021_01_09` | experiment date                                         |
| time                 | `numeric` | 0            | time / n generations, important for fitness calculation |
| group                | `numeric` | 1            | group number for sample                                 |
| reference_group      | `numeric` | 1            | group number of reference for comparison                |
| baseMean             | `numeric` | `NA`         | DESeq2 average number of reads for sgRNA                |
| log2FoldChange       | `numeric` | 0            | DESeq2 log2 FC for sgRNA                                |
| lfcSE                | `numeric` | 0            | DESeq2 log2 FC error for sgRNA                          |
| stat                 | `numeric` | `NA`         | DESeq2 t statistic for sgRNA                            |
| pvalue               | `numeric` | 1            | DESeq2 p-value for sgRNA                                |
| padj                 | `numeric` | 1            | DESeq2 adjusted p-value for sgRNA                       |
| fitness              | `numeric` | 2.020183     | fitness for sgRNA                                       |
| sgRNA_index          | `numeric` | 4            | relative position of sgRNA                              |
| sgRNA_correlation    | `numeric` | 0.6247412    | correlation of sgRNA with others                        |
| sgRNA_efficiency     | `numeric` | 0.9893041    | relative repression efficiency of sgRNA                 |
| wmean_log2FoldChange | `numeric` | 0            | weighted mean log2 FC for gene                          |
| sd_log2FoldChange    | `numeric` | 0            | standard deviation of log2 FC for gene                  |
| wmean_fitness        | `numeric` | 1.777574     | weighted mean fitness for gene                          |
| sd_fitness           | `numeric` | 0.9558989    | standard dev of fitness for gene                        |
| p_fitness            | `numeric` | 0.001        | p-value from Wilcoxon rank sum test (Null: fitness ~ 0) |
| p_fitness_adj        | `numeric` | 0.0001       | p-value from Wilcoxon test, Benjamini-Hochberg adjusted |
| comb_score           | `numeric` | 0.0001       | combined score (`-log10(p-value) * abs(wmean_fitness)`  |

### R markdown report

<details markdown="1">
<summary>Output files</summary>

- `fitness_report/`
  - `counts_summary.nb.html`: HTML report with information about all samples and read counts.
  - `fitness_summary.nb.html`: HTML report with information about fitness scores for all conditions.

</details>

Custom R markdown templates are used to render two HTML reports with information about all samples, their number of mapped reads, barcodes, genes, fitness scores, and other information. If fitness calculation is omitted through the option `gene_fitness = false` in the call to `nextflow run ...`, only the first report is written. The same applies for having less than two time points in the `samplesheet.csv` table.

### MultiQC

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

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
