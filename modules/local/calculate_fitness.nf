process FITNESS {
    tag "$samplesheet"
    label "process_high"

    conda (params.enable_conda ? "conda-forge::r-base=4.0 conda-forge::r-tidyverse bioconda::bioconductor-deseq2=1.28.0 bioconda::bioconductor-biocparallel bioconda::bioconductor-limma" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-31ad840d814d356e5f98030a4ee308a16db64ec5:0e852a1e4063fdcbe3f254ac2c7469747a60e361-0' :
        'rocker/tidyverse' }"

    input:
    path samplesheet
    path counts
    val normalization
    val gene_fitness
    val gene_sep

    output:
    path 'all_counts.tsv', emit: allcounts
    path 'result.Rdata', emit: rdata
    path 'result.tsv', emit: tsv
    path 'versions.yml', emit: versions

    script: // This script is bundled with the pipeline, in nf-core/crispriscreen/bin/
    def args = task.ext.args ?: ''
    def gene_controls = (params.gene_controls == null) ? '' : params.gene_controls

    """
    calculate_fitness.R \
        "${samplesheet}" \
        "${counts}" \
        "${normalization}" \
        "${gene_fitness}" \
        "${gene_sep}" \
        "${gene_controls}" \
        $task.cpus \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        bioconductor-limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """
}
