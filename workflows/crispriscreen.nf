/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowCrispriscreen.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Input fasta file for sgRNA library not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// MODULE: Local non nf-core modules
//
include { PREPARE_LIBRARY } from '../modules/local/prepare_library'
include { FITNESS } from '../modules/local/calculate_fitness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE                  } from '../modules/nf-core/trimgalore/main'
include { SEQTK_SAMPLE                } from '../modules/nf-core/seqtk/sample/main'
include { BOWTIE2_BUILD               } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN               } from '../modules/nf-core/bowtie2/align/main'
include { SUBREAD_FEATURECOUNTS       } from '../modules/nf-core/subread/featurecounts/main'
include { RMARKDOWNNOTEBOOK           } from '../modules/nf-core/rmarkdownnotebook/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow CRISPRISCREEN {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Convert *.fasta input file to pseudo genome database (*.saf)
    //
    PREPARE_LIBRARY (
        ch_fasta
    )
    ch_versions = ch_versions.mix(PREPARE_LIBRARY.out.versions)

    //
    // MODULE: Run Seqtk/sample (optional)
    //
    ch_reads = Channel.empty()
    if (params.subsampling) {
        SEQTK_SAMPLE (
            INPUT_CHECK.out.reads, params.sample_size
        )
        ch_reads = SEQTK_SAMPLE.out.reads
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    } else {
        ch_reads = INPUT_CHECK.out.reads
    }

    //
    // MODULE: Run Trim galore to cut adapters and filter by quality
    //
    ch_trimmedreads = Channel.empty()
    if (params.skip_trimming) {
        ch_trimmedreads = ch_reads
    } else {
        TRIMGALORE (
            ch_reads
        )
        ch_trimmedreads = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_trimmedreads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Bowtie2  - build genome database index from fasta input
    //
    BOWTIE2_BUILD (
        ch_fasta
    )
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Bowtie2  - align (filtered) reads to reference
    //
    BOWTIE2_ALIGN (
        ch_trimmedreads,
        BOWTIE2_BUILD.out.index,
        params.save_unaligned,
        params.sort_bam
    )
    ch_versions = ch_versions.mix(BOWTIE2_ALIGN.out.versions)

    //
    // MODULE: Subread/featureCounts - generates statistics about read counts per gene
    //
    ch_bowtiebam = Channel.empty()
    ch_bowtiebam = BOWTIE2_ALIGN.out.bam
    ch_bowtiebam
        .combine(PREPARE_LIBRARY.out.annotation)
        .set { ch_bowtiebam }

    SUBREAD_FEATURECOUNTS (
        ch_bowtiebam
    )
    ch_versions = ch_versions.mix(SUBREAD_FEATURECOUNTS.out.versions.first())

    //
    // MODULE: Calculate guide RNA and gene fitness score from read counts using DESeq2
    //
    ch_featurecounts = Channel.empty()
    ch_featurecounts = SUBREAD_FEATURECOUNTS.out.counts
    ch_featurecounts
        .map { it[1] }
        .collect()
        .set { ch_featurecounts }

    FITNESS (
        ch_input, ch_featurecounts, params.normalization,
        params.gene_fitness, params.gene_sep
    )
    ch_versions = ch_versions.mix(FITNESS.out.versions)

    //
    // MODULE: R markdown rendering the final fitness reports
    //
    ch_rmdtemplates = Channel.of(
        [ [ id:'counts_summary' ], file("$projectDir/bin/counts_summary.Rmd") ],
        [ [ id:'fitness_summary' ], file("$projectDir/bin/fitness_summary.Rmd") ]
    )
    ch_fitness = FITNESS.out.allcounts
    ch_fitness = ch_fitness.mix(FITNESS.out.rdata)
    ch_fitness = ch_fitness.collect()

    RMARKDOWNNOTEBOOK (
        ch_rmdtemplates, [], ch_fitness
    )
    ch_versions = ch_versions.mix(RMARKDOWNNOTEBOOK.out.versions)

    //
    // MODULE: Dump Software Versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowCrispriscreen.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))

    ch_multiqc_conflogo = Channel.of(
        [ ch_multiqc_config, file("https://multiqc.info/logos/MultiQC_logo.png") ]
    )
    MULTIQC (
        ch_multiqc_files.collect(), ch_multiqc_conflogo
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
