/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// // Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// // Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.fasta) { ch_fasta = file(params.fasta) } else { exit 1, 'Input fasta file for sgRNA library not specified!' }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Local non nf-core modules
//
include { PREPARE_LIBRARY } from '../modules/local/preparation/prepare_library'
include { PREPARE_COUNTS } from '../modules/local/preparation/prepare_counts'
include { FITNESS } from '../modules/local/fitness/calculate_fitness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { TRIMGALORE             } from '../modules/nf-core/trimgalore/main'
include { CUTADAPT               } from '../modules/nf-core/cutadapt/main'
include { SEQTK_SAMPLE           } from '../modules/nf-core/seqtk/sample/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN          } from '../modules/nf-core/bowtie2/align/main'
include { SUBREAD_FEATURECOUNTS  } from '../modules/nf-core/subread/featurecounts/main'
include { MAGECK_MLE             } from '../modules/nf-core/mageck/mle/main'
include { RMARKDOWNNOTEBOOK      } from '../modules/nf-core/rmarkdownnotebook/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_crispriscreen_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CRISPRISCREEN {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

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
            ch_samplesheet, params.sample_size
        )
        ch_reads = SEQTK_SAMPLE.out.reads
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
    } else {
        ch_reads = ch_samplesheet
    }

    //
    // MODULE: Run Trim galore to cut _generic_ (e.g. illumina) adapters and filter by quality
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
    // MODULE: Run cutadapt to cut _specific adapters_ such as primer sequences
    //
    ch_cutreads = Channel.empty()
    if (params.skip_trimming) {
        ch_cutreads = ch_trimmedreads
    } else {
        CUTADAPT (
            ch_trimmedreads
        )
        ch_cutreads = CUTADAPT.out.reads
        ch_versions = ch_versions.mix(CUTADAPT.out.versions)
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_cutreads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Bowtie2 - build genome database index from fasta input
    //
    ch_fasta_bowtie = [ [ id:'fasta' ], ch_fasta ]
    BOWTIE2_BUILD (
        ch_fasta_bowtie
    )
    ch_bowtie2_index = BOWTIE2_BUILD.out.index
    ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)

    //
    // MODULE: Bowtie2 - align (filtered) reads to reference
    //
    BOWTIE2_ALIGN (
        ch_cutreads,
        ch_bowtie2_index,
        ch_fasta_bowtie,
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

    ch_featurecounts = Channel.empty()
    ch_featurecounts = SUBREAD_FEATURECOUNTS.out.counts
    ch_featurecounts
        .map { it[1] }
        .collect()
        .set { ch_featurecounts }

    //
    // MODULE: Calculate gene fitness from read counts using Mageck
    //
    if (params.design_matrix) {
        ch_input_design = file(params.design_matrix, checkIfExists: true)
    } else {
        ch_input_design = file("|")
    }

    PREPARE_COUNTS (
        ch_input, ch_featurecounts, params.gene_sep, ch_input_design
    )

    ch_all_counts = PREPARE_COUNTS.out.all_counts
        .map { [ [ id: "all_counts" ], it ] }

    if (params.run_mageck & params.gene_fitness) {
        MAGECK_MLE (
            ch_all_counts,
            PREPARE_COUNTS.out.design
        )
        ch_versions = ch_versions.mix(MAGECK_MLE.out.versions)
    }

    //
    // MODULE: Calculate gene fitness from read counts using DESeq2
    //
    FITNESS (
        ch_input, PREPARE_COUNTS.out.all_counts, params.normalization,
        params.gene_fitness, params.gene_sep
    )
    ch_versions = ch_versions.mix(FITNESS.out.versions)

    //
    // MODULE: R markdown rendering the final fitness reports
    //
    rmd_templ_counts = [ [ id:'counts_summary' ], file("$projectDir/bin/counts_summary.Rmd") ]
    rmd_templ_fitness = [ [ id:'fitness_summary' ], file("$projectDir/bin/fitness_summary.Rmd") ]
    if (params.gene_fitness) {
        ch_rmdtemplates = Channel.of(
            rmd_templ_counts, rmd_templ_fitness
        )
    } else [
        ch_rmdtemplates = Channel.of(
            rmd_templ_counts
        )
    ]

    ch_fitness = FITNESS.out.rdata
    ch_fitness = ch_fitness.mix(PREPARE_COUNTS.out.all_counts)
    ch_fitness = ch_fitness.collect()

    RMARKDOWNNOTEBOOK (
        ch_rmdtemplates, [ test: 'somevalue'], ch_fitness
    )
    ch_versions = ch_versions.mix(RMARKDOWNNOTEBOOK.out.versions)

    //
    // MODULE: Dump Software Versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))
    ch_multiqc_files                      = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(TRIMGALORE.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(TRIMGALORE.out.log.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files                      = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.collect{it[1]}.ifEmpty([]))

    ch_multiqc_conflogo = Channel.of(
        [ ch_multiqc_config, file("https://multiqc.info/logos/MultiQC_logo.png") ]
    )
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
