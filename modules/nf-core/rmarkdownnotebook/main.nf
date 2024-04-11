include { dump_params_yml; indent_code_block } from "./parametrize"

process RMARKDOWNNOTEBOOK {
    tag "$meta.id"
    label 'process_low'

    // Override of default container with one that has latest R, yaml, rmakrdown, and tidyverse packages.
    // Fix changes by running "nf-core modules patch"
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://ghcr.io/MPUSP/nf-core-crispriscreen-fitness:1.0.0' :
        'oras://ghcr.io/MPUSP/nf-core-crispriscreen-fitness:1.0.0' }"

    input:
    tuple val(meta), path(notebook)
    val parameters
    path input_files

    output:
    tuple val(meta), path("*.html")              , emit: report
    tuple val(meta), path("*.parameterised.Rmd") , emit: parameterised_notebook, optional: true
    tuple val(meta), path ("artifacts/*")        , emit: artifacts, optional: true
    tuple val(meta), path ("session_info.log")   , emit: session_info
    path  "versions.yml"                         , emit: versions
    path  "*.png"                                , emit: png, optional: true
    path  "*.svg"                                , emit: svg, optional: true
    path  "*.pdf"                                , emit: pdf, optional: true
    path  "*_table.csv"                          , emit: csv, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def parametrize = (task.ext.parametrize == null) ?  true : task.ext.parametrize
    def implicit_params = (task.ext.implicit_params == null) ? true : task.ext.implicit_params
    def meta_params = (task.ext.meta_params == null) ? true : task.ext.meta_params

    // Dump parameters to yaml file.
    // Using a yaml file over using the CLI params because
    //  * no issue with escaping
    //  * allows to pass nested maps instead of just single values
    def params_cmd = ""
    def render_cmd = ""
    if (parametrize) {
        nb_params = [:]
        if (implicit_params) {
            nb_params["cpus"] = task.cpus
            nb_params["artifact_dir"] = "artifacts"
            nb_params["input_dir"] = "./"
        }
        if (meta_params) {
            nb_params["meta"] = meta
        }
        nb_params += parameters
        params_cmd = dump_params_yml(nb_params)
        render_cmd = """\
            params = yaml::read_yaml('.params.yml')

            # Instead of rendering with params, produce a version of the R
            # markdown with param definitions set, so the notebook itself can
            # be reused
            rmd_content <- readLines('${prefix}.Rmd')

            # Extract YAML content between the first two '---'
            start_idx <- which(rmd_content == "---")[1]
            end_idx <- which(rmd_content == "---")[2]
            rmd_yaml_content <- paste(rmd_content[(start_idx+1):(end_idx-1)], collapse = "\\n")
            rmd_params <- yaml::yaml.load(rmd_yaml_content)

            # Override the params
            rmd_params[['params']] <- modifyList(rmd_params[['params']], params)

            # Recursive function to add 'value' to list elements, except for top-level
            add_value_recursively <- function(lst, is_top_level = FALSE) {
                if (!is.list(lst)) {
                    return(lst)
                }

                lst <- lapply(lst, add_value_recursively)
                if (!is_top_level) {
                    lst <- list(value = lst)
                }
                return(lst)
            }

            # Reformat nested lists under 'params' to have a 'value' key recursively
            rmd_params[['params']] <- add_value_recursively(rmd_params[['params']], is_top_level = TRUE)

            # Convert back to YAML string
            updated_yaml_content <- as.character(yaml::as.yaml(rmd_params))

            # Remove the old YAML content
            rmd_content <- rmd_content[-((start_idx+1):(end_idx-1))]

            # Insert the updated YAML content at the right position
            rmd_content <- append(rmd_content, values = unlist(strsplit(updated_yaml_content, split = "\\n")), after = start_idx)

            writeLines(rmd_content, '${prefix}.parameterised.Rmd')

            # Render based on the updated file
            rmarkdown::render('${prefix}.parameterised.Rmd', output_file='${prefix}.html', envir = new.env())
        """
    } else {
        render_cmd = "rmarkdown::render('${prefix}.Rmd', output_file='${prefix}.html')"
    }

    """
    # Dump .params.yml heredoc (section will be empty if parametrization is disabled)
    ${indent_code_block(params_cmd, 4)}

    # Create output directory
    mkdir artifacts

    # Set parallelism for BLAS/MKL etc. to avoid over-booking of resources
    export MKL_NUM_THREADS="$task.cpus"
    export OPENBLAS_NUM_THREADS="$task.cpus"
    export OMP_NUM_THREADS="$task.cpus"

    # Work around  https://github.com/rstudio/rmarkdown/issues/1508
    # If the symbolic link is not replaced by a physical file
    # output- and temporary files will be written to the original directory.
    mv "${notebook}" "${notebook}.orig"
    cp -L "${notebook}.orig" "${prefix}.Rmd"

    # Render notebook
    Rscript - <<EOF
        ${indent_code_block(render_cmd, 8)}
        writeLines(capture.output(sessionInfo()), "session_info.log")
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rmarkdown: \$(Rscript -e "cat(paste(packageVersion('rmarkdown'), collapse='.'))")
    END_VERSIONS
    """
}
