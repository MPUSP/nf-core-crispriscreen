process FITNESS {
    tag "$fitness"
    label "process_high"

    conda (params.enable_conda ? "conda-forge::r-base=4.0 conda-forge::r-tidyverse bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel bioconda::bioconductor-limma" : null)

    input:
    path samplesheet
    val counts
    val normalization
    val gene_fitness
    val gene_sep

    output:
    path 'all_counts.tsv', emit: allcounts
    path 'result.Rdata', emit: rdata
    path 'result.tsv', emit: tsv
    path 'versions.yml', emit: versions

    script: // This script is bundled with the pipeline, in nf-core/crispriscreen/bin/
    """
    calculate_fitness.R \\
        $samplesheet \\
        $task.cpus \\
        "${counts}" \\
        "${normalization}" \\
        "${gene_fitness}" \\
        "${gene_sep}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        bioconductor-limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """
}
