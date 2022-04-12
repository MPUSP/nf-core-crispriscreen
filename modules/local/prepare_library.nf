process PREPARE_LIBRARY {
    tag "$prepare_libary"
    label "process_low"

    conda (params.enable_conda ? "conda-forge::r-base=4.0" : null)

    input:
    path fasta

    output:
    path '*.saf'       , emit: annotation
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/crispriscreen/bin/
    """
    prepare_library.R \\
        $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
    END_VERSIONS
    """
}
