process FITNESS {
    tag "$samplesheet"
    label "process_high"

    conda (params.enable_conda ? "conda-forge::r-base=4.2.2 conda-forge::r-tidyverse=1.3.2 bioconda::bioconductor-deseq2=1.38.0 bioconda::bioconductor-biocparallel=1.32.0 bioconda::bioconductor-limma=3.54.0" : null)

    input:
    path samplesheet
    path counts
    val normalization
    val gene_fitness
    val gene_sep

    output:
    path 'result.Rdata', emit: rdata, optional:true
    path 'result.tsv', emit: tsv, optional:true
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
        r-tidyverse: \$(Rscript -e "library(tidyverse); cat(as.character(packageVersion('tidyverse')))")
        bioconductor-deseq2: \$(Rscript -e "library(DESeq2); cat(as.character(packageVersion('DESeq2')))")
        bioconductor-biocparallel: \$(Rscript -e "library(biocparallel); cat(as.character(packageVersion('biocparallel')))")
        bioconductor-limma: \$(Rscript -e "library(limma); cat(as.character(packageVersion('limma')))")
    END_VERSIONS
    """
}
